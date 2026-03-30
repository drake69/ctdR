# ctdR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/drake69/ctdR/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/drake69/ctdR/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/drake69/ctdR/actions/workflows/test-coverage.yaml/badge.svg?branch=main)](https://github.com/drake69/ctdR/actions/workflows/test-coverage.yaml)
[![Codecov](https://codecov.io/gh/drake69/ctdR/branch/main/graph/badge.svg)](https://codecov.io/gh/drake69/ctdR)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![GitHub issues](https://img.shields.io/github/issues/drake69/ctdR)](https://github.com/drake69/ctdR/issues)
[![GitHub stars](https://img.shields.io/github/stars/drake69/ctdR)](https://github.com/drake69/ctdR/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/drake69/ctdR)](https://github.com/drake69/ctdR/commits/main)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![Bioconductor dependencies](https://img.shields.io/badge/Bioconductor-dependencies-green.svg)](https://www.bioconductor.org/)
<!-- badges: end -->

**ctdR** is an R package that identifies chemicals significantly associated with a set of genes using data from the [Comparative Toxicogenomics Database (CTD)](https://ctdbase.org).

## Features

- **Over-Representation Analysis (ORA)** — tests whether the overlap between your gene list and each chemical's target genes is larger than expected by chance (via [`clusterProfiler`](https://bioconductor.org/packages/clusterProfiler/))
- **Gene Set Enrichment Analysis (GSEA)** — detects chemicals whose targets cluster toward the extremes of a ranked gene list (via [`fgsea`](https://bioconductor.org/packages/fgsea/))
- **Automatic caching** — CTD data is parsed once and cached locally for fast repeated analyses
- **Human-only filtering** — automatically restricts interactions to *Homo sapiens* (OrganismID 9606)

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

### From GitHub

```r
# install.packages("devtools")
devtools::install_github("drake69/ctdR")
```

### Bioconductor dependencies

ctdR depends on several Bioconductor packages. If they are not installed automatically, run:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("fgsea", "org.Hs.eg.db", "clusterProfiler", "DOSE", "AnnotationDbi"))
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

### Step 3 — Run enrichment analysis

```r
# Prepare your gene list as a data frame
genes <- data.frame(
  entrez_ids = c("7124", "3569", "7157", "672", "1956"),
  pvalue     = c(0.001, 0.003, 0.01, 0.02, 0.05)
)

# Over-Representation Analysis (default)
ora_results <- enrichment_CTD(genes, method = "ORA")
head(ora_results)

# Gene Set Enrichment Analysis
gsea_results <- enrichment_CTD(genes, method = "GSEA")
head(gsea_results)
```

## Input Format

The `entrez_ids` parameter must be a data frame with at least two columns:

| Column | Description |
|---|---|
| `entrez_ids` | Character or numeric Entrez gene IDs |
| *(second column)* | A numeric value per gene (e.g. p-value, log fold-change). Used for ranking in GSEA; ignored in ORA. |

## Output

### ORA results

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

### GSEA results

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
| `foldEnrichment` | \|ES\| / mean(ES) |
| `Enriched_GENE` | Comma-separated enriched gene symbols |

## Dependencies

### CRAN

- [readr](https://cran.r-project.org/package=readr) — fast CSV reading
- [rappdirs](https://cran.r-project.org/package=rappdirs) — cross-platform cache directory
- [plyr](https://cran.r-project.org/package=plyr) — data manipulation

### Bioconductor

- [fgsea](https://bioconductor.org/packages/fgsea/) — fast GSEA implementation
- [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) — ORA enrichment
- [DOSE](https://bioconductor.org/packages/DOSE/) — ratio parsing utilities
- [AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi/) — gene ID mapping
- [org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/) — human gene annotation

## Contributing

Contributions are welcome! Please open an [issue](https://github.com/drake69/ctdR/issues) or submit a pull request.

## Citation

If you use ctdR in your research, please cite:

```
Corsaro L (2024). ctdR: Enrichment Analysis of Chemical-Gene Interactions
from the Comparative Toxicogenomics Database. R package version 0.1.2.
https://github.com/drake69/ctdR
```

Additionally, please cite CTD as required by their [Terms of Service](https://ctdbase.org/about/legal.jsp):

> Davis AP, Wiegers TC, Johnson RJ, Sciaky D, Wiegers J, Mattingly CJ.
> Comparative Toxicogenomics Database (CTD): update 2023.
> *Nucleic Acids Research*. 2023;51(D1):D1257-D1262.

## License

Apache License 2.0 - see [LICENSE](LICENSE) for details.
