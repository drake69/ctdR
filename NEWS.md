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
