# Changes in version 0.99.4

## Breaking changes

These changes are made now, while the package is pre-1.0 and not yet
accepted into Bioconductor, so that the public column schema is stable
before any external code depends on it.

### Input column name (ORA / GSEA)

* The input data frame must now provide an `EntrezID` column (was
  `entrez_ids`). The numeric value column can still be named freely.

### Unified output schema (ORA, GSEA, CAMERA)

All three data-frame-returning methods now share the same leading
columns, in this order:

```
ChemicalID, ChemicalName, Method, PValue, PValueAdjusted, ...
```

The new `Method` column carries the method label (`"ORA"`, `"GSEA"`,
`"CAMERA"`), making cross-method `rbind` / `dplyr::bind_rows`
straightforward.

Rows are sorted by `PValueAdjusted` ascending in all three methods
(previously ORA results were unsorted).

### Column renames (full table)

| Method | Old column | New column |
|---|---|---|
| GSEA | `pval` | `PValue` |
| GSEA | `ES` | `EnrichmentScore` |
| GSEA | `NES` | `NormalizedEnrichmentScore` |
| GSEA | `size` | `GeneSetSize` |
| GSEA | `leadingEdge` | `LeadingEdge` |
| GSEA | `Enriched_GENE` | `EnrichedGenes` |
| ORA  | `pvalue` | `PValue` |
| ORA  | `padj` | `PValueAdjusted` |
| ORA  | `BgRatio` | `BackgroundRatio` |
| ORA  | `qvalue` | `QValue` |
| ORA  | `geneID` | `EnrichedGenes` |
| ORA  | `foldEnrichment` | `FoldEnrichment` |
| ORA  | `ID` | `ChemicalID` |
| ORA  | `Description` | *(dropped — was always `== ID`)* |
| CAMERA | `PValue` | (kept as) `PValue` |
| CAMERA | `NGenes` | `GeneSetSize` |
| CAMERA | `FDR` | *(dropped — `PValueAdjusted` recomputed per `pAdjustMethod`)* |
| All  | (new) | `Method` |

Cross-method semantic alignment:
- `GeneSetSize` replaces both GSEA's `size` and CAMERA's `NGenes`
- `EnrichedGenes` replaces both ORA's `geneID` and GSEA's `Enriched_GENE`
- `PValue` / `PValueAdjusted` are spelled the same across all methods

## Internal refactor

* `enrichment_CTD()` shrank from 101 lines to ~30 by delegating
  argument validation to `.validate_enrichment_args()` and adding
  `.run_ora()` / `.run_gsea()` runners mirroring the existing
  `.run_camera()` / `.run_gsva()` shape.
* `.format_enrichment_result()` is the **single source of truth** for
  the engine -> canonical-schema mapping: each runner passes a
  `rename = c(old = "New")` and `drop = c(...)` and the formatter
  handles padj recomputation, metadata merge, column ordering, sort
  and row-name reset.
* `gsea()` and `ora()` engines now return their underlying tool's
  native column casing (`pval` / `pvalue`, `ES`, `NES`, …); the only
  semantic rename they apply themselves is lifting the primary-key
  column (`pathway` for fgsea, `ID` for clusterProfiler) to
  `ChemicalID`, so internal callers never see the misleading generic
  name.
* `gsea()` signature dropped the unused `chemicals` and
  `pAdjustMethod` arguments. Its `entrez_ids` parameter was renamed
  to `gene_table` (still internal-only).
* `.run_camera()` shrank to 48 lines (was 57) and `.run_ora()` /
  `.run_gsea()` are now under 30 lines each. BiocCheck's "function
  length > 50" NOTE is satisfied for the enrichment-table pipeline.

# Changes in version 0.99.2

## New features

* `enrichment_CTD()` now supports four enrichment methods through a unified
  interface, selectable via the `method` argument:
  - `"ORA"` (default) — Over-Representation Analysis via
    `clusterProfiler::enricher()` (unchanged).
  - `"GSEA"` — rank-based Gene Set Enrichment Analysis via `fgsea::fgsea()`
    (unchanged).
  - `"CAMERA"` — competitive gene-set test accounting for inter-gene
    correlation, via `limma::camera()`. Input: a numeric expression matrix
    (genes x samples) plus a design matrix and a contrast.
  - `"GSVA"` — per-sample Gene Set Variation Analysis via `GSVA::gsva()`,
    returning a chemical x sample score matrix.
* New first argument `x` (polymorphic): a data frame for ORA/GSEA, a numeric
  matrix for CAMERA/GSVA. Auto-detects identifier type (Entrez vs HGNC SYMBOL)
  from `rownames(x)`; an explicit `id_type` override is also accepted.
* `plot_CTD()` now dispatches on input class and method-specific columns:
  bar/dot plots of fold enrichment for ORA/GSEA, bar/dot plots of
  `-log10(padj)` coloured by direction of enrichment for CAMERA, and a
  sample-level heatmap of the top-variance chemicals for GSVA.

## Deprecation

* The first argument of `enrichment_CTD()` was renamed `entrez_ids` -> `x`.
  Calls using the old name still work but emit a deprecation warning and will
  be removed in a future release.

## Dependencies

* New `Imports`: `limma`, `GSVA`, `stats`.

# Changes in version 0.99.1

## Bioconductor reviewer feedback

* Removed `renv` from the package. Dependencies are managed via
  `DESCRIPTION` and installed by `r-lib/actions/setup-r-dependencies`
  (pak-based) in CI.

## Bug fixes

* `import_CTD()` now stores `ChemicalName_GeneSymbols$gene` as a
  character vector instead of a factor, fixing
  "universe must be a character vector" in `clusterProfiler::enricher()`.
* Inlined `parse_ratio()` as a private helper because `DOSE::parse_ratio`
  is no longer exported. `DOSE` is dropped from `Imports`.
* Aligned `man/gsea.Rd` parameter name with `R/gsea.R`
  (`ChemicalName_GeneEntrezIds`), removing a codoc-mismatch WARNING.
* Declared `plot_CTD()` ggplot2 NSE column references in
  `utils::globalVariables()`, removing "no visible binding" NOTEs.

## Infrastructure

* Pinned `trufflesecurity/trufflehog` to a concrete version; the bare
  `@v3` tag does not exist in the action repo.
* Made `oysteR` audit fail-soft when OSS Index credentials are missing
  (still fails the build on real vulnerabilities).
* Removed the `dependency-review` job: GitHub's Dependency Graph does
  not support R `DESCRIPTION` files.

# Changes in version 0.99.0

## Improvements

* Bumped version to 0.99.0 for Bioconductor submission.
* Fixed R CMD check to pass with 0 errors, 0 warnings, 0 notes.
* Fixed CI workflows for macOS, Ubuntu, and Windows.
* Added Codecov integration for coverage reporting.
* Added GitHub Pages site with usage examples.
* Fixed broken README badges.

# Changes in version 0.1.2

## New features

* Added `import_CTD()` function to import and cache CTD chemical-gene
  interaction data from user-downloaded files.
* Added Over-Representation Analysis (ORA) method via `clusterProfiler::enricher`.
* `enrichment_CTD()` now supports two methods: `"ORA"` (default) and `"GSEA"`.
* Clear error message when CTD data has not been imported, with instructions
  on how to download and import the required file.

## Improvements

* Separated data import from analysis — users call `import_CTD()` once, then
  `enrichment_CTD()` as many times as needed.
* Added fold enrichment calculation to both ORA and GSEA results.
* `gsea()` now receives `chemicals` as an explicit parameter instead of
  relying on the parent environment.
* Removed deprecated `nperm` parameter from `fgsea::fgsea` call.
* Removed leftover `browser()` calls.

## Documentation

* Added comprehensive roxygen documentation for all exported and internal
  functions, including parameter descriptions, return value tables, and
  examples.
* Added package-level help page (`?ctdR`) with quick start guide.
* Added data licensing disclaimer throughout documentation and DESCRIPTION.
* Added README.md with badges, installation instructions, and usage examples.
* Added vignette with complete workflow.

## Testing

* Added testthat test suite with 100% line coverage.
* Tests cover: import validation, missing data errors, ORA enrichment,
  GSEA enrichment, NA handling, and column structure.

## Infrastructure

* Added GitHub Actions CI for macOS, Ubuntu, and Windows.
* Added test coverage workflow with Codecov integration.

# Changes in version 0.1.1

* Added GSEA analysis via `fgsea::fgsea`.
* Initial caching mechanism using `rappdirs`.

# Changes in version 0.1.0

* Initial package skeleton.
* Basic gene-chemical interaction parsing from CTD CSV files.
