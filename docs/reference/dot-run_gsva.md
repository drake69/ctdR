# Per-sample gene-set scoring with GSVA

Internal function that runs Gene Set Variation Analysis
([`gsva`](https://rdrr.io/pkg/GSVA/man/gsva.html)) on an expression
matrix using the cached CTD chemical gene sets. Unlike ORA/GSEA/CAMERA —
which return a single p-value per chemical (group-level inference) —
GSVA produces a matrix of **per-sample enrichment scores**: one row per
chemical, one column per sample. These scores are suitable for
downstream clustering, association tests against phenotypes, survival
analysis, or heatmap visualization.

## Usage

``` r
.run_gsva(expr, id_type = NULL, cache_dir, ...)
```

## Arguments

- expr:

  Numeric expression matrix (genes x samples). `rownames(expr)` must be
  Entrez IDs or HGNC symbols matching the cached CTD gene sets.

- id_type:

  Either `"entrez"`, `"symbol"`, or `NULL` for auto-detection from
  `rownames(expr)`.

- cache_dir:

  Directory holding the cached CTD `.rda` files.

- ...:

  Forwarded to
  [`gsvaParam`](https://rdrr.io/pkg/GSVA/man/gsvaParam-class.html) (e.g.
  `kcdf`, `minSize`, `maxSize`, `tau`, `maxDiff`).

## Value

A numeric matrix of GSVA enrichment scores with CTD chemical IDs in rows
and samples in columns.
