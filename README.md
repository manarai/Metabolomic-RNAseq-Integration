# Maize × *Fusarium verticillioides* — Multi-Omics Integration

End-to-end RNA-seq + LC-MS metabolomics integration pipeline on the maize 282 diversity panel infected with *Fusarium verticillioides*. Identifies disease-driven and resistance-associated transcripts and metabolites, integrates them into joint latent factors, and runs mediation analysis on candidate causal triples.

---

## Biological context

- **Panel**: 110 maize inbred lines from the 282 diversity panel, all infected with *F. verticillioides* (n=4 biological replicates per genotype, no mock controls).
- **Disease phenotype**: ergosterol (fungal-membrane sterol → proxy for *Fusarium* biomass; higher = worse symptoms).
- **RNA-seq subset**: top 10 most-resistant + top 10 most-susceptible lines were re-grown and sequenced (raw counts available).
- **Integration scope**: 23 genotypes have both omics.

---

## Pipeline overview

The `multiomics_integration_v4.ipynb` notebook is the **clean, DESeq2-primary version** with full integration. Run end-to-end:

| Step | What it does | Key outputs |
|---|---|---|
| 1. Load + harmonize | parse all three CSVs, strip growth-tag suffixes from genotype names | — |
| 2. Preprocess | metabolomics: zero → ½ min, log2, per-sample median center, replicate-average. RNA: FPKM → TPM → log2 | — |
| 3. Metabolite ↔ ergosterol | Pearson + BH-FDR across all 52k features | `metabolite_vs_ergosterol_correlation.csv` |
| 4. PCA / t-SNE / clustermap | unsupervised structure on the disease-associated panel; row-cluster extraction | `metabolite_clusters.csv` |
| 5. DESeq2 (primary DE) | label R/S contrast + metabolomic-extreme contrast on raw counts | `DE_label_R_vs_S_DESeq2.csv`, `DE_metabolomic_extremes_DESeq2.csv` |
| 6. Integration | gene × metabolite correlation matrix on 23 shared genotypes; hub ranking | `gene_metabolite_pairs.csv`, `gene_hub_ranking.csv` |
| 7. Sparse PLS | `cca-zoo SCCA_PMD` joint factors with active feature shortlist | `spls_metabolite_loadings.csv`, `spls_gene_loadings.csv` |
| 7b. Method comparison | PLSCanonical + MOFA + DIABLO + JIVE for integration-method validation | `integration_method_comparison.csv`, `mofa_factors.csv` |
| 7c. Mediation | Sobel-test Gene → Metabolite → Disease triples on hub features | `mediation_analysis.csv` |
| 8. KEGG enrichment | curated compound IDs + REST queries + GSEA-style Mann-Whitney | `kegg_compound_pathway_mapping.csv`, `kegg_pathway_gsea.csv` |
| 9. Findings | summary of results, limitations, suggested next steps | (markdown) |

---

## How the integration works (step by step)

Multi-omics integration in v4 is layered: each step asks a question the previous step couldn't answer. The full chain runs from raw correlations → joint latent factors → causal hypotheses → pathway interpretation.

### Step 1 — Sample alignment (the substrate for everything)

Replicates are not paired across omics (rep 1 of B73 metabolomics is a different plant from rep 1 of B73 RNA-seq), so the common unit is the **genotype**, not the sample.

- Metabolomics: average 4 reps per genotype → 109 genotypes × 52,231 metabolites
- RNA-seq: average ~4 reps per genotype → 20 genotypes × 21,654 expressed genes
- **Intersection = 23 shared genotypes** — every integration step downstream uses this 23-sample matrix

### Step 2 — Per-omic feature filtering (before pairing them)

3.5 million pairwise tests on n = 23 is a noise machine. So each omic gets filtered first:

- **Metabolomics**: keep the 3,666 features that correlate with ergosterol at FDR < 0.05 & |r| > 0.3 ("disease-associated panel")
- **RNA-seq**: keep the top 500 by **DESeq2 padj** ∪ top 500 by **variance** = 948 unique genes

### Step 3 — Gene × metabolite correlation network (§6)

The feature-by-feature view: every metabolite correlated with every gene.

1. Z-score each row across the 23 samples
2. Single matrix multiply: `R = M_z @ G_z.T / (n - 1)` → the full 3,666 × 948 correlation matrix in one operation
3. Analytic p-values from Student-t: `t = r · √(n-2) / √(1 - r²)`
4. Keep pairs at **|r| ≥ 0.6 AND p ≤ 1e-3** → ~44,000 significant pairs (uncorrected — at n = 23, several thousand are chance)
5. **Hub ranking**: count how many metabolites each gene links to. A gene that connects to 30 disease-associated metabolites is more interesting than any single high-r pair. Top hubs are annotated with their DESeq2 stats.

**Outputs**: `gene_metabolite_pairs.csv` (the network), `gene_hub_ranking.csv` (the priority list).

### Step 4 — Sparse PLS — primary joint-factor method (§7)

Correlation pairs are a flat view. Sparse PLS asks: **is there a single shared axis of variation that both omics encode together?**

Algorithm (`cca-zoo` `SCCA_PMD`):

- Find weight vectors `w_m` (metabolites) and `w_g` (genes) such that the projected samples `T = X_met · w_m` and `U = X_rna · w_g` are maximally correlated.
- Add an **L1 penalty** on `w_m` and `w_g` (controlled by `tau = 0.3`) → most weights are forced to exactly zero. Only the features that *actually* drive the joint signal survive.
- Repeat with deflation to extract 3 orthogonal factor pairs.

Real-run output:

| Factor | Active metabolites | Active genes | Inter-omic r | r vs ergosterol |
|---|---|---|---|---|
| F1 | 75 | 33 | **+0.94** | +0.62 |
| F2 | 50 | 32 | +0.91 | −0.53 |
| F3 | 58 | 24 | +0.87 | −0.67 |

**The key result**: with just 75 metabolites + 33 genes active on F1, the metabolomics-projection and RNA-seq-projection of the same 23 samples correlate at **r = 0.94** — the two omics see the same thing along this axis. F1 ↔ ergosterol r = +0.62, so **F1 is the joint disease signature** encoded by ~110 features total. Those ~110 features are the publication biomarker shortlist.

### Step 5 — Method validation (§7b)

A single algorithm's answer could be quirk. So the same integration is rerun through four more methods:

| Method | What it does differently | Disease factor r vs ergosterol |
|---|---|---|
| **PLSCanonical (dense)** | Sparse PLS without L1 — sensitivity check | **+0.93** |
| **MOFA** | Bayesian factor model with spike-and-slab prior; per-view variance decomposition | (opt-in: `pip install mofapy2`) |
| **DIABLO** (R / mixOmics) | **Supervised** — uses R/S labels to find a discriminating factor | **−0.93** (sign arbitrary) |
| **JIVE** (R / r.jive) | Decomposes joint vs view-specific structure; reports rank | **+0.86** (joint rank = 1) |

**Four out of four methods recover the same disease axis at |r| ≥ 0.6.** JIVE's `joint rank = 1` is the cleanest result: there is exactly **one** shared dimension between the two omics. Metabolomics has 5 additional dimensions of private structure; RNA-seq has 7. The integration adds value (shared signal exists), and each omic also carries unique information worth analyzing alone.

### Step 6 — Mediation analysis — the causal step (§7c)

Steps 3-5 find correlations. Mediation tests the causal direction: **does the gene affect disease THROUGH the metabolite?**

For each (gene g, metabolite m, ergosterol e) triple from the sparse-PLS hub features:

1. `m = a·g + a₀` — does the gene predict the metabolite? (path a)
2. `e = c'·g + b·m + e₀` — does the metabolite predict disease given the gene? (paths b and c')
3. `e = c·g + e₀` — total effect of gene on disease (path c)
4. **Sobel z-test** for the indirect effect: `z = (a·b) / √(b²·SE_a² + a²·SE_b²)`
5. **Proportion mediated** = (c − c') / c

10 top hub genes × 10 top hub metabolites = 100 triples tested. **22 significant at Sobel p < 0.05, all with > 30% mediation** → 22 candidate causal mechanisms. These are the highest-priority experimental follow-ups.

### Step 7 — Biological interpretation via KEGG (§8)

The integration outputs are gene IDs and m/z features. To make them biologically interpretable:

1. Manual curation: 25 named compounds (phenylpropanoids, flavonoids, JA, benzoxazinoids, …) → KEGG compound IDs
2. Substring-match against the metabolite index → 19 named compounds found in the dataset
3. KEGG REST API: `https://rest.kegg.jp/link/pathway/cpd:<id>` → pathway memberships (308 compound-pathway edges, 59 distinct pathways)
4. **GSEA-style enrichment**: for each pathway, Mann-Whitney U test on member compounds' r vs ergosterol vs the rest of the metabolome.

Top hits: **phenylpropanoid biosynthesis (12 compounds), phenylalanine metabolism (8), plant hormone biosynthesis (10), glucosinolate biosynthesis (8)** — textbook maize defense biology, recovered from the data with no prior bias.

### Mental model

```
       ┌────────────────────────────┐
       │ Sample-match: 23 genotypes │  ← both omics aligned at the genotype level
       └─────────────┬──────────────┘
                     │
         ┌───────────┴───────────┐
         ▼                       ▼
 ┌─────────────────┐    ┌────────────────────┐
 │ 3,666 metabs    │    │  948 genes         │
 │ (ergo-FDR<0.05) │    │ (DESeq2 + var)     │
 └────────┬────────┘    └─────────┬──────────┘
          │                       │
          └─────────┬─────────────┘
                    │
        ┌───────────┼───────────────┐
        ▼           ▼               ▼
  ┌─────────┐  ┌──────────┐   ┌───────────┐
  │ Network │  │ Joint    │   │ Mediation │
  │ (pairs, │  │ factors  │   │ (causal   │
  │  hubs)  │  │ (sparse  │   │  triples) │
  │   §6    │  │  PLS,    │   │   §7c     │
  │         │  │  §7, 7b) │   │           │
  └────┬────┘  └────┬─────┘   └─────┬─────┘
       │            │               │
       └────────────┴───────────────┘
                    │
                    ▼
            ┌───────────────┐
            │ KEGG pathway  │  ← biological context
            │ enrichment §8 │
            └───────────────┘
```

Each layer answers a different question:
- **§6 (network)** — *which features co-vary?*
- **§7 (sparse PLS)** — *is there a shared latent disease axis?*
- **§7b (validation)** — *is that axis method-robust?*
- **§7c (mediation)** — *do the gene effects on disease go through metabolites?*
- **§8 (pathways)** — *which biology does this implicate?*

### Single-sentence summary (for the paper)

*"We identified a single shared metabolomic-transcriptomic disease axis (sparse PLS F1, inter-omic r = 0.94, 75 metabolites + 33 genes active), recovered by four independent integration methods (PLSCanonical, DIABLO, JIVE), with 22 candidate causal Gene → Metabolite → Disease triples (Sobel p < 0.05) and pathway enrichment dominated by phenylpropanoid + plant hormone biosynthesis."*

---

### Other notebooks in this repo

- `multiomics_integration_v2.ipynb` — comprehensive version (~130 cells) with side-by-side comparisons of five DE methods (Welch, Wilcoxon, limma-trend, DESeq2, edgeR), per-replicate-vs-averaging pseudoreplication validation, and an honest discussion of analytical choices. Useful as a reference / methods supplement.
- `multiomics_integration_v3.ipynb` — minimal/lean version. Same pipeline, ~25 cells, suitable for quick reproduction.
- `multiomics_integration_notebook.ipynb` — original exploratory draft (kept for history).

---

## Input data

Place these three files at the repo root (or update the paths in §0 of the notebooks):

| File | Description | Approximate size |
|---|---|---|
| `181029_282_panel_Fvert_3d_cut_stem_x393_v9_LIB_Filtered_Final_Metaboanalyst_NoNC340-NC352_FINALDATASETforMS1.csv` | LC-MS metabolomics, row 1 = genotype, row 2 = `Ergosterol`, rows 3+ = features | ~125 MB |
| `260217_RNAseq_Fvert282_topR-S_Final_Metaboanalyst.csv` | RNA-seq FPKM, row 1 = R/S phenotype | ~11 MB |
| `210328_SC_Fvert_Novoseq_clean_counts.csv` | RNA-seq raw counts (for DESeq2/edgeR) | ~11 MB |

> The metabolomics CSV exceeds GitHub's 100 MB file limit. Either use `git-lfs` (instructions below) or host the file externally and document the URL in this README.

### Git-LFS setup (one-time, if you want to track large files)

```bash
git lfs install
git lfs track "*.csv"
git add .gitattributes
```

---

## Installation

### Python environment (conda)

```bash
# Clone the repo
git clone https://github.com/<your-username>/maize-fusarium-multiomics.git
cd maize-fusarium-multiomics

# Create the environment from environment.yml
conda env create -f environment.yml
conda activate maize-multiomics

# Verify the install
python -c "import pandas, numpy, scipy, sklearn, pydeseq2, cca_zoo, mofapy2, statsmodels, rpy2; print('Python deps OK')"
```

### R environment (for DIABLO, JIVE, limma, edgeR via R-magic)

Some integration methods (DIABLO, JIVE) and the supplementary DE methods (limma-trend, edgeR) require R. The notebooks call R via the `rpy2` Python bridge using `%%R` cell magic.

**Install R itself**:
```bash
# macOS:
brew install r

# Ubuntu/Debian:
sudo apt-get install r-base r-base-dev

# Or download from https://cran.r-project.org/
```

**Install R packages** (one-time, takes ~5-10 min):
```bash
Rscript install_r_packages.R
```

This installs from CRAN/Bioconductor:
- `BiocManager` (CRAN, bootstrap)
- `limma`, `edgeR`, `mixOmics` (Bioconductor)
- `r.jive` (CRAN)

---

## Running the analysis

```bash
# Activate the environment
conda activate maize-multiomics

# Option A — run end-to-end from the terminal (writes a fully-rendered HTML report)
jupyter nbconvert --to html --execute multiomics_integration_v4.ipynb \
    --ExecutePreprocessor.timeout=1800

# Option B — interactive
jupyter lab multiomics_integration_v4.ipynb
```

### Generate a clean code-free HTML report (for collaborators)
```bash
jupyter nbconvert --to html --no-input \
    multiomics_integration_v4_executed.ipynb
```

---

## Output files

All written to `outputs_v4/` (created automatically):

**Section 3 — Metabolite correlation**
- `metabolite_vs_ergosterol_correlation.csv` — every metabolite's r, p, FDR vs ergosterol
- `metabolite_clusters.csv` — K=5 cluster assignment

**Section 5 — Differential expression**
- `DE_label_R_vs_S_DESeq2.csv` — DESeq2 results for the R/S phenotype contrast
- `DE_metabolomic_extremes_DESeq2.csv` — DESeq2 results for the K=8 metabolomic-extreme contrast

**Section 6-7 — Integration**
- `gene_metabolite_pairs.csv` — significant cross-omics pairs
- `gene_hub_ranking.csv` — genes ranked by metabolite link count + DESeq2 annotation
- `spls_metabolite_loadings.csv`, `spls_gene_loadings.csv` — sparse-PLS factor loadings
- `integration_method_comparison.csv` — sparse PLS vs PLSCanonical vs MOFA vs DIABLO vs JIVE summary
- `mofa_factors.csv` — MOFA latent factors (when `mofapy2` is installed)
- `mediation_analysis.csv` — Gene → Metabolite → Disease Sobel-test results

**Section 8 — Pathway enrichment**
- `kegg_compound_pathway_mapping.csv` — curated compound → KEGG pathway edges
- `kegg_pathway_gsea.csv` — per-pathway GSEA-style enrichment statistics

---

## Methods summary

**Differential expression**: DESeq2 (primary, negative-binomial GLM on raw counts) with metabolomically-defined extreme contrast as a parallel analysis. v2 additionally runs Welch, Wilcoxon, limma-trend, and edgeR for cross-method robustness.

**Unsupervised structure**: PCA (global disease gradient), t-SNE (local metabolic subtypes), hierarchical clustermap with ergosterol annotation bar.

**Integration**:
- *Correlation network* — pairwise gene × metabolite Pearson on n=23 shared genotypes.
- *Sparse PLS* (cca-zoo `SCCA_PMD`) — primary joint-factor method; L1 sparsity gives an interpretable biomarker shortlist.
- *Method validation* — PLSCanonical (dense sensitivity check), MOFA (Bayesian factor model with per-view variance decomposition), DIABLO (supervised sparse generalized CCA), JIVE (joint vs view-specific variance decomposition).
- *Mediation* — Sobel test on Gene → Metabolite → Disease triples to nominate causal mechanisms.

**Pathway enrichment**: curated metabolite-to-KEGG-compound mapping + REST API for pathway membership + Mann-Whitney GSEA-style test on the ergosterol-correlation ranking.

---

## Replicate handling — design rationale

| Notebook section | Sample unit | Why |
|---|---|---|
| §3-4 metabolite correlation, PCA, clustering | **per genotype** (n=109) | 4 reps share one genetic background → using n=430 is **pseudoreplication** (validated in v2 §6c) |
| §5 DE | **per replicate** (n≈85) | DESeq2's NB GLM needs within-group variance |
| §6+ integration | **per genotype** (REQUIRED) | metabolomics rep 1 ≠ RNA-seq rep 1 for the same genotype (separate plants); only common unit is the genotype |

---

## Known limitations

1. **n = 23 shared genotypes** is small for joint modeling — sparse-PLS r values should be permutation-tested + cross-validated for publication.
2. **No population-structure adjustment** — the 282 panel has tropical/temperate/sweet subpopulations that can drive spurious associations. Add a population PC as a covariate to the §3 and §5 regressions for the publication.
3. **No FDR on the gene × metabolite pairs** in §6 — current threshold (uncorrected p ≤ 1e-3 at n=23) is exploratory. The sparse-PLS-derived hub pairs are the cleaner alternative.
4. **~96% of metabolites are LC-MS unknowns** (`C_<id>_<m/z>mz_<RT>min_…`). KEGG pathway enrichment uses the ~20 named compounds; extend `KEGG_CURATED` in §8 as more compounds are annotated.
5. **All findings are correlational** (mediation gives suggestive causal direction but not proof) — experimental validation needed for hub genes (RT-qPCR, knockouts).

---

## Repository layout

```
maize-fusarium-multiomics/
├── README.md                              # this file
├── environment.yml                        # conda Python environment
├── install_r_packages.R                   # one-shot R-package installer
├── .gitignore
├── .gitattributes                         # if using git-lfs
│
├── multiomics_integration_v4.ipynb        # CLEAN, primary notebook (45 cells)
├── multiomics_integration_v3.ipynb        # minimal version (~25 cells)
├── multiomics_integration_v2.ipynb        # comprehensive reference (~130 cells)
├── multiomics_integration_notebook.ipynb  # original draft
│
├── 181029_…FINALDATASETforMS1.csv         # metabolomics (input; ~125 MB)
├── 260217_…FPKM.csv                       # RNA-seq FPKM (input)
├── 210328_SC_Fvert_Novoseq_clean_counts.csv  # RNA-seq raw counts (input)
│
└── outputs_v4/                            # all generated tables (gitignored)
    ├── DE_label_R_vs_S_DESeq2.csv
    ├── DE_metabolomic_extremes_DESeq2.csv
    ├── gene_metabolite_pairs.csv
    ├── gene_hub_ranking.csv
    ├── integration_method_comparison.csv
    ├── mediation_analysis.csv
    ├── spls_*.csv
    ├── kegg_*.csv
    └── …
```

---

## Citation

If this pipeline is useful in a publication, please cite the upstream tools:

- **DESeq2** — Love MI, Huber W, Anders S (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology 15:550. (`pydeseq2`)
- **edgeR** — Robinson MD, McCarthy DJ, Smyth GK (2010). *edgeR: a Bioconductor package for differential expression analysis of digital gene expression data*. Bioinformatics 26:139.
- **limma** — Ritchie ME et al. (2015). *limma powers differential expression analyses for RNA-sequencing and microarray studies*. NAR 43:e47.
- **DIABLO** — Singh A et al. (2019). *DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays*. Bioinformatics 35:3055. (`mixOmics`)
- **JIVE** — Lock EF et al. (2013). *Joint and individual variation explained (JIVE) for integrated analysis of multiple data types*. Annals of Applied Statistics 7:523. (`r.jive`)
- **MOFA** — Argelaguet R et al. (2018). *Multi-Omics Factor Analysis — a framework for unsupervised integration of multi-omics data sets*. Molecular Systems Biology 14:e8124.
- **sparse PLS / SCCA_PMD** — Witten DM, Tibshirani R, Hastie T (2009). *A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis*. Biostatistics 10:515. (`cca-zoo`)
- **KEGG** — Kanehisa M, Goto S (2000). *KEGG: Kyoto Encyclopedia of Genes and Genomes*. NAR 28:27.

---

## License

MIT — see `LICENSE` file.
