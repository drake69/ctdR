# Gene Set Enrichment Analysis (GSEA)

Internal function that performs Gene Set Enrichment Analysis using
[`fgsea`](https://rdrr.io/pkg/fgsea/man/fgsea.html). Detects chemicals
whose known target genes cluster toward the extremes of a ranked gene
list.

## Usage

``` r
gsea(ChemicalName_GeneEntrezIds, entrez_ids, chemicals, pAdjustMethod = "BH")
```

## Arguments

- ChemicalName_GeneEntrezIds:

  A named list where each element is a character vector of Entrez gene
  IDs associated with a CTD chemical. Names are CTD chemical IDs.

- entrez_ids:

  A data frame with at least two columns:

  `entrez_ids`

  :   Entrez gene IDs (character).

  (second column)

  :   Numeric values used for ranking (e.g. p-values). Converted
      internally to `-log10(value)` for the ranking statistic.

- chemicals:

  A data frame with columns `ChemicalID` and `ChemicalName`, used to
  annotate results with human-readable names.

- pAdjustMethod:

  Character. Method for multiple testing correction (default `"BH"`).
  Passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) to
  recalculate the `padj` column returned by `fgsea`.

## Value

A data frame with columns: `ChemicalID`, `pval`, `padj`, `ES`, `NES`,
`size`, `leadingEdge`, `ChemicalName`, `foldEnrichment`, and
`Enriched_GENE`. Results are sorted by `padj` in ascending order.
