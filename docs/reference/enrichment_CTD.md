# Chemical Enrichment Analysis Using CTD

Identifies chemicals whose known gene targets are significantly enriched
in a user-supplied gene list or expression matrix, using data from the
Comparative Toxicogenomics Database (CTD).

Four methods are available:

- **ORA**:

  Over-Representation Analysis (default). Tests whether the overlap
  between your gene list and each chemical's target genes is larger than
  expected by chance. Uses
  [`enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).
  Input: data frame of Entrez IDs (column `entrez_ids`) with an optional
  numeric value column.

- **GSEA**:

  Gene Set Enrichment Analysis. Uses a ranked gene list (ranked by the
  numeric column in the input, e.g. p-values or fold changes) to detect
  chemicals whose targets cluster toward the top or bottom of the
  ranking. Uses [`fgsea`](https://rdrr.io/pkg/fgsea/man/fgsea.html).

- **CAMERA**:

  Competitive gene-set test accounting for inter-gene correlation. Uses
  [`camera`](https://rdrr.io/pkg/limma/man/camera.html). Input: a
  numeric expression matrix (genes x samples) plus a design matrix and a
  contrast.

- **GSVA**:

  Gene Set Variation Analysis (sample-level scoring). Returns a matrix
  of per-sample enrichment scores for each chemical, suitable for
  downstream clustering or association testing. Uses
  [`gsva`](https://rdrr.io/pkg/GSVA/man/gsva.html). Input: a numeric
  expression matrix (genes x samples).

## Usage

``` r
enrichment_CTD(
  x,
  method = c("ORA", "GSEA", "CAMERA", "GSVA"),
  design = NULL,
  contrast = NULL,
  id_type = NULL,
  pAdjustMethod = "BH",
  ...
)
```

## Arguments

- x:

  The input. Its expected type depends on `method`:

  - For `"ORA"` and `"GSEA"`: a data frame with at least two columns,
    `entrez_ids` (character or numeric Entrez gene IDs) and a numeric
    value column (e.g. p-value, log fold-change). The numeric column is
    used for ranking in GSEA and ignored in ORA.

  - For `"CAMERA"` and `"GSVA"`: a numeric expression matrix with genes
    in rows and samples in columns. `rownames(x)` must be either Entrez
    IDs or HGNC SYMBOLs.

- method:

  Character. Enrichment method: `"ORA"` (default), `"GSEA"`, `"CAMERA"`,
  or `"GSVA"`.

- design:

  Design matrix (required when `method = "CAMERA"`).

- contrast:

  Contrast specification for
  [`camera`](https://rdrr.io/pkg/limma/man/camera.html) (column number,
  column name, or numeric vector). Required when `method = "CAMERA"`.

- id_type:

  Either `"entrez"`, `"symbol"`, or `NULL` (default) for auto-detection
  from `rownames(x)`. Only used when `method` is `"CAMERA"` or `"GSVA"`.

- pAdjustMethod:

  Character. Method for multiple testing correction: one of `"BH"`
  (Benjamini-Hochberg, default), `"bonferroni"`, `"fdr"` (alias for BH),
  or `"none"`. Not used for `method = "GSVA"` (which returns scores
  rather than p-values).

- ...:

  Additional arguments forwarded to the underlying engine:
  [`enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
  for ORA, [`fgsea`](https://rdrr.io/pkg/fgsea/man/fgsea.html) for GSEA,
  [`camera`](https://rdrr.io/pkg/limma/man/camera.html) for CAMERA,
  [`gsva`](https://rdrr.io/pkg/GSVA/man/gsva.html) for GSVA.

## Value

- For `"ORA"`, `"GSEA"`, and `"CAMERA"`: a data frame of enrichment
  results sorted by adjusted p-value (`padj`).

- For `"GSVA"`: a numeric matrix of enrichment scores with chemicals
  (CTD chemical IDs) in rows and samples in columns.

## Details

Before calling this function you must import the CTD data once with
[`import_CTD`](https://drake69.github.io/ctdR/reference/import_CTD.md).
If the cached data is not found, the function stops with an informative
error message.

## Data Licensing Disclaimer

This package does **not** bundle or redistribute any CTD data. The
Comparative Toxicogenomics Database is maintained by NC State University
and its data are subject to specific licensing terms. Users must
download the data directly from <https://ctdbase.org> and comply with
the CTD Terms of Service (<https://ctdbase.org/about/legal.jsp>).

## See also

[`import_CTD`](https://drake69.github.io/ctdR/reference/import_CTD.md)
to import and cache the CTD data;
[`plot_CTD`](https://drake69.github.io/ctdR/reference/plot_CTD.md) to
visualize results.

## Examples

``` r
# Import the bundled sample data first:
sample_file <- system.file(
    "extdata", "CTD_chem_gene_ixns_sample.csv",
    package = "ctdR"
)
import_CTD(sample_file)
#> Reading CTD chemical-gene interactions from: /Library/Frameworks/R.framework/Versions/4.6/Resources/library/ctdR/extdata/CTD_chem_gene_ixns_sample.csv
#> Filtered to 86 human interactions
#> Mapping genes for 10 chemicals (this may take a while)...
#> 
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> CTD data cached successfully in: ~/Library/Caches/ctdR

# ORA / GSEA: prepare a gene list with Entrez IDs and a numeric value
genes <- data.frame(
    entrez_ids = c("7124", "3569", "7157", "672", "1956"),
    pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
)
ora_results <- enrichment_CTD(genes, method = "ORA")
#> 'select()' returned 1:1 mapping between keys and columns
gsea_results <- enrichment_CTD(genes, method = "GSEA")
#> Warning: All values in the stats vector are greater than zero and scoreType is "std", maybe you should switch to scoreType = "pos".

# CAMERA: expression matrix + design + contrast
set.seed(1)
expr <- matrix(rnorm(20 * 6), nrow = 20,
    dimnames = list(c("7124","3569","7157","672","1956",
        "1017","1019","207","208","595",
        "894","983","1029","1869","5599",
        "5290","4609","6347","3458","2099"),
        paste0("S", 1:6)))
grp <- factor(rep(c("ctrl","treat"), each = 3))
d <- model.matrix(~ grp)
# \donttest{
camera_results <- enrichment_CTD(expr, method = "CAMERA",
    design = d, contrast = 2)

# GSVA: per-sample enrichment scores
gsva_scores <- enrichment_CTD(expr, method = "GSVA")
#> ℹ GSVA version 2.6.1
#> ℹ Searching for rows with constant values
#> ℹ Calculating GSVA ranks
#> ℹ kcdf='auto' (default)
#> ℹ GSVA dense (classical) algorithm
#> ℹ Row-wise ECDF estimation with Gaussian kernels
#> ℹ Calculating row ECDFs
#> ℹ Calculating column ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Calculating GSVA scores for 10 gene sets
#> ✔ Calculations finished
# }
```
