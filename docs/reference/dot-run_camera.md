# Competitive gene-set enrichment with CAMERA

Internal function that runs the CAMERA competitive gene-set test
([`camera`](https://rdrr.io/pkg/limma/man/camera.html)) of an expression
matrix against the cached CTD chemical gene sets. Unlike ORA/GSEA —
which take a single (ranked) gene list — CAMERA tests, for each
chemical, whether its target genes show a stronger differential signal
than the rest of the transcriptome under a user-supplied design and
contrast, while accounting for inter-gene correlation.

## Usage

``` r
.run_camera(
  expr,
  design,
  contrast,
  id_type = NULL,
  pAdjustMethod = "BH",
  chemicals_meta,
  cache_dir,
  ...
)
```

## Arguments

- expr:

  Numeric expression matrix (genes x samples). `rownames(expr)` must be
  Entrez IDs or HGNC symbols matching the cached CTD gene sets.

- design:

  Design matrix produced e.g. by
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html).

- contrast:

  Either a column number or column name of `design`, or a numeric
  contrast vector.

- id_type:

  Either `"entrez"`, `"symbol"`, or `NULL` for auto-detection from
  `rownames(expr)`.

- pAdjustMethod:

  Method passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)
  to compute the `padj` column.

- chemicals_meta:

  Data frame with columns `ChemicalID` and `ChemicalName`, used to
  annotate results.

- cache_dir:

  Directory holding the cached CTD `.rda` files.

- ...:

  Forwarded to [`camera`](https://rdrr.io/pkg/limma/man/camera.html)
  (e.g. `inter.gene.cor`, `use.ranks`, `allow.neg.cor`, `trend.var`).

## Value

A data frame with columns `ChemicalID`, `ChemicalName`, `NGenes`,
`Direction`, `Correlation` (if reported), `pvalue`, `padj`, sorted by
`padj` ascending.
