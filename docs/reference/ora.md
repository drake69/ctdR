# Over-Representation Analysis (ORA)

Internal function that performs Over-Representation Analysis using
[`enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).
Tests whether the input gene list is significantly enriched for genes
associated with each CTD chemical.

## Usage

``` r
ora(ChemicalName_GeneSymbols, entrez_ids, pAdjustMethod = "BH")
```

## Arguments

- ChemicalName_GeneSymbols:

  A data frame with two columns (`term`, `gene`) mapping CTD chemical
  IDs to HGNC gene symbols. This serves as the TERM2GENE input for
  [`clusterProfiler::enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).

- entrez_ids:

  Character vector of HGNC gene symbols to test for enrichment.

- pAdjustMethod:

  Character. Method for multiple testing correction (default `"BH"`).
  Passed to
  [`enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).

## Value

A data frame with columns: `ID`, `GeneRatio`, `BgRatio`, `pvalue`,
`padj`, `qvalue`, `geneID`, `Count`, and `foldEnrichment`. Returns an
empty data frame with the same structure if no enrichment is found.
