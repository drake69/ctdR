# ctdR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/drake69/ctdR/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/drake69/ctdR/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/drake69/ctdR/actions/workflows/test-coverage.yaml/badge.svg?branch=main)](https://github.com/drake69/ctdR/actions/workflows/test-coverage.yaml)
[![Codecov](https://codecov.io/gh/drake69/ctdR/branch/main/graph/badge.svg)](https://codecov.io/gh/drake69/ctdR)
[![security-scan](https://github.com/drake69/ctdR/actions/workflows/security.yaml/badge.svg?branch=main)](https://github.com/drake69/ctdR/actions/workflows/security.yaml)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19344200.svg)](https://doi.org/10.5281/zenodo.19344200)
[![GitHub issues](https://img.shields.io/github/issues/drake69/ctdR)](https://github.com/drake69/ctdR/issues)
[![GitHub stars](https://img.shields.io/github/stars/drake69/ctdR)](https://github.com/drake69/ctdR/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/drake69/ctdR)](https://github.com/drake69/ctdR/commits/main)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.5-blue.svg)](https://www.r-project.org/)
[![Bioconductor dependencies](https://img.shields.io/badge/Bioconductor-dependencies-green.svg)](https://www.bioconductor.org/)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://drake69.github.io/ctdR/)
<!-- badges: end -->

**ctdR** is an R package that identifies chemicals significantly associated with a set of genes using data from the [Comparative Toxicogenomics Database (CTD)](https://ctdbase.org).

> 📖 **Documentation, tutorials & function reference**: **[drake69.github.io/ctdR](https://drake69.github.io/ctdR/)** — full package website with articles and examples.

> ⭐ **Find ctdR useful?** [**Star it on GitHub**](https://github.com/drake69/ctdR/stargazers) — it takes one second and helps other researchers discover the package.

> 🚧 **Bioconductor status**: *not yet accepted, or under review* — submission [#4232](https://github.com/Bioconductor/Contributions/issues/4232). *Fingers crossed* 🤞. For now, install from GitHub (instructions below).

## Features

Four enrichment methods through a unified `enrichment_CTD()` interface, selected via the `method` argument:

- **ORA**: Over-Representation Analysis (hypergeometric test). Input: gene list. Backend: [`stats::phyper`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html), computed directly.
- **GSEA** — Gene Set Enrichment Analysis (rank-based, permutational). Input: ranked gene list. Backend: [`fgsea::fgsea`](https://bioconductor.org/packages/fgsea/).
- **CAMERA**: competitive gene-set test with inter-gene correlation correction. Input: expression matrix or `SummarizedExperiment`, plus design + contrast. Backend: [`limma::camera`](https://bioconductor.org/packages/limma/).
- **GSVA**: Gene Set Variation Analysis (per-sample scoring). Input: expression matrix or `SummarizedExperiment`. Output: chemical × sample scores, returned in whichever of the two you supplied. Backend: [`GSVA::gsva`](https://bioconductor.org/packages/GSVA/).

Additional features:

- **Automatic caching** — CTD data is parsed once and cached locally for fast repeated analyses
- **Human-only filtering** — automatically restricts interactions to *Homo sapiens* (OrganismID 9606)
- **Auto-detected gene identifiers**: Entrez vs HGNC SYMBOL detected from `rownames()` of the input matrix or `SummarizedExperiment` (CAMERA / GSVA); override via `id_type`
- **Visualization** — `plot_CTD()` auto-dispatches to bar/dot plot (ORA/GSEA/CAMERA) or heatmap (GSVA)

## Data Licensing Disclaimer

> **This package does NOT bundle, redistribute, or embed any data from the Comparative Toxicogenomics Database.**
>
> CTD data are created and maintained by [NC State University](https://ncsu.edu/) and are subject to specific licensing terms and conditions. Users are solely responsible for:
>
> 1. Downloading the required data files directly from <https://ctdbase.org>
> 2. Reading and accepting the [CTD Terms of Service](https://ctdbase.org/about/legal.jsp) before using the data
> 3. Complying with all applicable CTD data licensing requirements, including proper citation of CTD in any publications or derived works
>
> By using ctdR with CTD data, you acknowledge that you have read, understood, and agreed to the CTD Terms of Service.

## Installation

> **Status**: ctdR is currently under review for inclusion in Bioconductor
> ([Bioconductor/Contributions#4232](https://github.com/Bioconductor/Contributions/issues/4232))
> — *not on Bioconductor yet, fingers crossed* 🤞. Until then, please
> install from GitHub.

### From GitHub (current installation method)

```r
# install.packages("devtools")
devtools::install_github("drake69/ctdR")
```

### Bioconductor dependencies

ctdR depends on several Bioconductor packages. They are pulled in
automatically by `install_github()`, but you can install them upfront with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("fgsea", "org.Hs.eg.db",
                       "AnnotationDbi", "limma", "GSVA"))
```

### From Bioconductor (once accepted)

Once ctdR is accepted into Bioconductor (currently under review — *fingers
crossed* 🤞), it will be installable directly via:

```r
BiocManager::install("ctdR")
```

## Quick Start

### Step 1 — Download CTD data (once)

Download **`CTD_chem_gene_ixns.csv.gz`** from the CTD website:

<https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz>

Decompress the file:

```bash
gunzip CTD_chem_gene_ixns.csv.gz
```

### Step 2 — Import into ctdR (once)

```r
library(ctdR)

import_CTD("~/Downloads/CTD_chem_gene_ixns.csv")
#> Reading CTD chemical-gene interactions from: ~/Downloads/CTD_chem_gene_ixns.csv
#> Filtered to 1234567 human interactions
#> Mapping genes for 12345 chemicals (this may take a while)...
#> CTD data cached successfully in: /Users/you/Library/Caches/ctdR
```

The data is now cached locally. You only need to do this once (or again when you download a newer CTD release).

`import_CTD()` also accepts a URL in place of a local path; remote files are downloaded and cached via `BiocFileCache`, with a one-time reminder pointing to the CTD data-licensing terms. The package assumes no default URL, so you always supply the source yourself.

### Step 3 — Run enrichment analysis

The first argument of `enrichment_CTD()` is polymorphic and named `x`:
a `data.frame` for ORA / GSEA, and for CAMERA / GSVA either a numeric
matrix or a [`SummarizedExperiment`](https://bioconductor.org/packages/SummarizedExperiment).

#### ORA / GSEA — gene-list paradigm

```r
# Prepare your gene list as a data frame (Entrez IDs + numeric column)
genes <- data.frame(
  entrez_ids = c("7124", "3569", "7157", "672", "1956"),
  pvalue     = c(0.001, 0.003, 0.01, 0.02, 0.05)
)

# Over-Representation Analysis (default)
ora_results <- enrichment_CTD(genes, method = "ORA")
head(ora_results)

# Gene Set Enrichment Analysis (uses the second column as ranking)
gsea_results <- enrichment_CTD(genes, method = "GSEA")
head(gsea_results)

# Customize multiple testing correction (default is "BH")
ora_bonf <- enrichment_CTD(genes, method = "ORA", pAdjustMethod = "bonferroni")
```

#### CAMERA — multi-sample paradigm with correlation correction

```r
# Suppose `se` is a SummarizedExperiment whose assay holds normalised
# expression (genes x samples) with Entrez IDs (or HGNC SYMBOLs) as
# rownames, and whose colData carries the experimental group.
# A plain numeric matrix works too; the SummarizedExperiment keeps the
# sample annotation attached to the data.
design <- model.matrix(~ se$group)

camera_results <- enrichment_CTD(
  se,
  method   = "CAMERA",
  design   = design,
  contrast = 2  # last column of `design` = treat vs ctrl
)
head(camera_results)
```

#### GSVA — per-sample scoring

```r
# SummarizedExperiment in, SummarizedExperiment out: the scores arrive
# with the sample annotation still attached. Pass a matrix instead and
# you get a plain matrix back.
gsva_scores <- enrichment_CTD(se, method = "GSVA")
dim(gsva_scores)
SummarizedExperiment::colData(gsva_scores)
```

### Step 4 — Visualize results

```r
# ORA / GSEA / CAMERA: bar or dot plot of top chemicals
plot_CTD(ora_results,    type = "bar")
plot_CTD(gsea_results,   type = "dot", n = 10)
plot_CTD(camera_results, type = "bar")  # colour encodes Direction

# GSVA: heatmap of the top-variance chemicals across samples
plot_CTD(gsva_scores)
```

### Utilities and interoperability

```r
# Read a CTD file into a DataFrame via the Bioconductor import() convention
library(BiocIO)
ctd <- import(CTDFile("~/Downloads/CTD_chem_gene_ixns.csv"))

# Retrieve processed tables from the cache
chem <- ctd_cache("chemicals")

# Export the CTD gene sets for a third-party engine (e.g. EnrichmentBrowser)
gene_sets <- as_genesets_CTD("entrez")
```

In `ctdR`, **CTD** always means the *Comparative Toxicogenomics Database*. It is unrelated to the CRAN package `CTD` ("Connect The Dots") and to the `ctd` (CellTypeDataset) object used by `EWCE`.

## End-to-end runnable example

For a production-shaped reference pipeline — full GSE311566 Female
PBMCs dataset (Dex vs DMSO), a-priori power analysis, declared
significance thresholds, `limma` differential expression, and all
four enrichment methods with BH-adjusted alpha cutoffs — see:

```
inst/scripts/example_gse311566_full_pipeline.R
```

Run it from the package root after installing ctdR and populating
the CTD cache (`import_CTD()`):

```bash
Rscript inst/scripts/example_gse311566_full_pipeline.R
```

Outputs (DE table, per-method significant chemicals, plots, power
report, `sessionInfo`) land in `./example_outputs/`. This script is
the example referenced by the companion paper; the vignette covers
the same workflow on a bundled subset (~34 KB, no alpha cutoffs)
for didactic purposes and BiocCheck-friendly build times.

## Input Format

The first argument `x` of `enrichment_CTD()` depends on the method:

### ORA / GSEA — data frame

| Column | Description |
|---|---|
| `entrez_ids` | Character or numeric Entrez gene IDs |
| *(second column)* | A numeric value per gene (e.g. p-value, log fold-change). Used for ranking in GSEA; ignored in ORA. |

### CAMERA / GSVA: numeric matrix or `SummarizedExperiment`

A `genes × samples` numeric matrix with `rownames(x)` set to either
Entrez IDs or HGNC SYMBOLs, or a
[`SummarizedExperiment`](https://bioconductor.org/packages/SummarizedExperiment)
wrapping one. Identifier type is auto-detected from `rownames(x)`;
override with the `id_type` argument if needed.

Prefer the `SummarizedExperiment` when you have sample annotation. It
keeps the assay and the annotation in one object, so subsetting or
reordering samples moves both together instead of leaving you to keep
a separate group vector in step with the columns. Use `assay` to pick
which assay to read when the object carries more than one; the default
is the first.

## Output

### ORA results (`data.frame`)

| Column | Description |
|---|---|
| `ChemicalID` | CTD chemical identifier (e.g. `"D000082"`) |
| `ChemicalName` | Human-readable chemical name |
| `GeneRatio` | Proportion of input genes in the chemical's set |
| `BgRatio` | Background ratio |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (BH method) |
| `foldEnrichment` | GeneRatio / BgRatio |
| `geneID` | Enriched gene symbols |
| `Count` | Number of overlapping genes |

### GSEA results (`data.frame`)

| Column | Description |
|---|---|
| `ChemicalID` | CTD chemical identifier |
| `ChemicalName` | Human-readable chemical name |
| `pval` | Raw p-value |
| `padj` | Adjusted p-value |
| `ES` | Enrichment score |
| `NES` | Normalized enrichment score |
| `size` | Size of the gene set |
| `leadingEdge` | Leading-edge gene subset |
| `foldEnrichment` | abs(ES) / mean(ES) |
| `Enriched_GENE` | Comma-separated enriched gene symbols |

### CAMERA results (`data.frame`)

| Column | Description |
|---|---|
| `ChemicalID` | CTD chemical identifier |
| `ChemicalName` | Human-readable chemical name |
| `NGenes` | Number of gene-set genes matched in `rownames(x)` |
| `Direction` | `"Up"` or `"Down"` — direction of mean effect |
| `Correlation` | Estimated (or pre-specified) inter-gene correlation |
| `pvalue` | Raw p-value from `limma::camera()` |
| `padj` | Adjusted p-value |

### GSVA results (`matrix` or `SummarizedExperiment`)

CTD chemical IDs in rows and samples in columns, containing per-sample
enrichment scores (typically in `[-1, 1]`). The container follows the
input: a matrix in returns a numeric matrix, a `SummarizedExperiment`
in returns a `SummarizedExperiment` whose assay holds the scores and
whose `colData` is carried over, so the annotation needed to interpret
the scores stays next to them.
Unlike the other methods, GSVA does not return p-values — the scores
are descriptive features for downstream analyses (clustering,
association with outcomes, heatmap visualization).

## Dependencies

Mirrors the `Imports` field of `DESCRIPTION`.

### CRAN

- [ggplot2](https://cran.r-project.org/package=ggplot2) — publication-quality plots
- [readr](https://cran.r-project.org/package=readr) — fast CSV reading

### Bioconductor

- [fgsea](https://bioconductor.org/packages/fgsea/) — fast GSEA implementation
- [limma](https://bioconductor.org/packages/limma/) — CAMERA backend (`limma::camera()`)
- [GSVA](https://bioconductor.org/packages/GSVA/) — per-sample gene-set scoring
- [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/): standard container accepted by CAMERA / GSVA
- [BiocIO](https://bioconductor.org/packages/BiocIO/): `import()` convention behind `CTDFile`
- [BiocFileCache](https://bioconductor.org/packages/BiocFileCache/): local cache of the processed CTD tables
- [S4Vectors](https://bioconductor.org/packages/S4Vectors/): `DataFrame` returned by `import()`
- [AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi/) — gene ID mapping
- [org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/) — human gene annotation

## Continuous Integration & Security

ctdR uses GitHub Actions for continuous integration and automated security auditing:

| Workflow | Purpose |
|---|---|
| **R-CMD-check** | Package build & check on macOS, Ubuntu, Windows |
| **test-coverage** | Code coverage via `covr` + Codecov |
| **security-scan** | Automated cybersecurity pipeline (see below) |

The **security-scan** workflow runs on every push/PR and weekly, and includes four jobs:

1. **Dependency vulnerability audit** — scans all installed R packages against the [Sonatype OSS Index](https://ossindex.sonatype.org/) via [`oysteR`](https://cran.r-project.org/package=oysteR), flagging packages with known CVEs.
2. **Static code analysis** — runs [`lintr`](https://cran.r-project.org/package=lintr) on the entire package source to detect code quality and potential security issues.
3. **Dependency review** (PR only) — uses GitHub's [dependency-review-action](https://github.com/actions/dependency-review-action) to flag new dependencies with high-severity vulnerabilities before merging.
4. **Secret & credential scan** — uses [TruffleHog](https://github.com/trufflesecurity/trufflehog) to detect accidentally committed secrets or API keys in the repository history.

## Contributing

Contributions are welcome! Please open an [issue](https://github.com/drake69/ctdR/issues) or submit a pull request.

## Citation

If you use ctdR in your research, please cite:

```
Corsaro L (2026). ctdR: Enrichment Analysis of Chemical-Gene Interactions
from the Comparative Toxicogenomics Database. R package version 0.99.2.
doi: 10.5281/zenodo.19344200. https://github.com/drake69/ctdR
```

Additionally, please cite CTD as required by their [Terms of Service](https://ctdbase.org/about/legal.jsp):

> Davis AP, Wiegers TC, Johnson RJ, Sciaky D, Wiegers J, Mattingly CJ.
> Comparative Toxicogenomics Database (CTD): update 2023.
> *Nucleic Acids Research*. 2023;51(D1):D1257-D1262.

## ⭐ If you find ctdR useful, please star the repo

A star helps other researchers discover the package and lets us see who's using it. Click the ⭐ at the [top of the GitHub page](https://github.com/drake69/ctdR/stargazers) — it takes one second and means a lot.

## License

Apache License 2.0 - see [LICENSE](LICENSE) for details.
