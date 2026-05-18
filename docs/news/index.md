# Changelog

## Changes in version 0.99.2

### New features

- [`enrichment_CTD()`](https://drake69.github.io/ctdR/reference/enrichment_CTD.md)
  now supports four enrichment methods through a unified interface,
  selectable via the `method` argument:
  - `"ORA"` (default) — Over-Representation Analysis via
    [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
    (unchanged).
  - `"GSEA"` — rank-based Gene Set Enrichment Analysis via
    [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html)
    (unchanged).
  - `"CAMERA"` — competitive gene-set test accounting for inter-gene
    correlation, via
    [`limma::camera()`](https://rdrr.io/pkg/limma/man/camera.html).
    Input: a numeric expression matrix (genes x samples) plus a design
    matrix and a contrast.
  - `"GSVA"` — per-sample Gene Set Variation Analysis via
    [`GSVA::gsva()`](https://rdrr.io/pkg/GSVA/man/gsva.html), returning
    a chemical x sample score matrix.
- New first argument `x` (polymorphic): a data frame for ORA/GSEA, a
  numeric matrix for CAMERA/GSVA. Auto-detects identifier type (Entrez
  vs HGNC SYMBOL) from `rownames(x)`; an explicit `id_type` override is
  also accepted.
- [`plot_CTD()`](https://drake69.github.io/ctdR/reference/plot_CTD.md)
  now dispatches on input class and method-specific columns: bar/dot
  plots of fold enrichment for ORA/GSEA, bar/dot plots of `-log10(padj)`
  coloured by direction of enrichment for CAMERA, and a sample-level
  heatmap of the top-variance chemicals for GSVA.

### Deprecation

- The first argument of
  [`enrichment_CTD()`](https://drake69.github.io/ctdR/reference/enrichment_CTD.md)
  was renamed `entrez_ids` -\> `x`. Calls using the old name still work
  but emit a deprecation warning and will be removed in a future
  release.

### Dependencies

- New `Imports`: `limma`, `GSVA`, `stats`.

## Changes in version 0.99.1

### Bioconductor reviewer feedback

- Removed `renv` from the package. Dependencies are managed via
  `DESCRIPTION` and installed by `r-lib/actions/setup-r-dependencies`
  (pak-based) in CI.

### Bug fixes

- [`import_CTD()`](https://drake69.github.io/ctdR/reference/import_CTD.md)
  now stores `ChemicalName_GeneSymbols$gene` as a character vector
  instead of a factor, fixing “universe must be a character vector” in
  [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).
- Inlined `parse_ratio()` as a private helper because
  `DOSE::parse_ratio` is no longer exported. `DOSE` is dropped from
  `Imports`.
- Aligned `man/gsea.Rd` parameter name with `R/gsea.R`
  (`ChemicalName_GeneEntrezIds`), removing a codoc-mismatch WARNING.
- Declared
  [`plot_CTD()`](https://drake69.github.io/ctdR/reference/plot_CTD.md)
  ggplot2 NSE column references in
  [`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html),
  removing “no visible binding” NOTEs.

### Infrastructure

- Pinned `trufflesecurity/trufflehog` to a concrete version; the bare
  `@v3` tag does not exist in the action repo.
- Made `oysteR` audit fail-soft when OSS Index credentials are missing
  (still fails the build on real vulnerabilities).
- Removed the `dependency-review` job: GitHub’s Dependency Graph does
  not support R `DESCRIPTION` files.

## Changes in version 0.99.0

### Improvements

- Bumped version to 0.99.0 for Bioconductor submission.
- Fixed R CMD check to pass with 0 errors, 0 warnings, 0 notes.
- Fixed CI workflows for macOS, Ubuntu, and Windows.
- Added Codecov integration for coverage reporting.
- Added GitHub Pages site with usage examples.
- Fixed broken README badges.

## Changes in version 0.1.2

### New features

- Added
  [`import_CTD()`](https://drake69.github.io/ctdR/reference/import_CTD.md)
  function to import and cache CTD chemical-gene interaction data from
  user-downloaded files.
- Added Over-Representation Analysis (ORA) method via
  [`clusterProfiler::enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).
- [`enrichment_CTD()`](https://drake69.github.io/ctdR/reference/enrichment_CTD.md)
  now supports two methods: `"ORA"` (default) and `"GSEA"`.
- Clear error message when CTD data has not been imported, with
  instructions on how to download and import the required file.

### Improvements

- Separated data import from analysis — users call
  [`import_CTD()`](https://drake69.github.io/ctdR/reference/import_CTD.md)
  once, then
  [`enrichment_CTD()`](https://drake69.github.io/ctdR/reference/enrichment_CTD.md)
  as many times as needed.
- Added fold enrichment calculation to both ORA and GSEA results.
- [`gsea()`](https://drake69.github.io/ctdR/reference/gsea.md) now
  receives `chemicals` as an explicit parameter instead of relying on
  the parent environment.
- Removed deprecated `nperm` parameter from
  [`fgsea::fgsea`](https://rdrr.io/pkg/fgsea/man/fgsea.html) call.
- Removed leftover [`browser()`](https://rdrr.io/r/base/browser.html)
  calls.

### Documentation

- Added comprehensive roxygen documentation for all exported and
  internal functions, including parameter descriptions, return value
  tables, and examples.
- Added package-level help page
  ([`?ctdR`](https://drake69.github.io/ctdR/reference/ctdR-package.md))
  with quick start guide.
- Added data licensing disclaimer throughout documentation and
  DESCRIPTION.
- Added README.md with badges, installation instructions, and usage
  examples.
- Added vignette with complete workflow.

### Testing

- Added testthat test suite with 100% line coverage.
- Tests cover: import validation, missing data errors, ORA enrichment,
  GSEA enrichment, NA handling, and column structure.

### Infrastructure

- Added GitHub Actions CI for macOS, Ubuntu, and Windows.
- Added test coverage workflow with Codecov integration.

## Changes in version 0.1.1

- Added GSEA analysis via
  [`fgsea::fgsea`](https://rdrr.io/pkg/fgsea/man/fgsea.html).
- Initial caching mechanism using `rappdirs`.

## Changes in version 0.1.0

- Initial package skeleton.
- Basic gene-chemical interaction parsing from CTD CSV files.
