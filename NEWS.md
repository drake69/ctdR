# Changes in version 0.99.9

## Significant user-visible changes

* ORA no longer goes through `clusterProfiler::enricher()`. The
  hypergeometric test is computed directly with `stats::phyper()`, and
  `clusterProfiler` has been removed from `Imports`. The p-values are
  unchanged: `scripts/ora_equivalence_check.R` in the revisions
  repository compares the two implementations across 24 configurations
  and finds them identical. The reason is cost, not correctness.
  `clusterProfiler` accounted for 59 of the package's 175 hard
  dependencies for this single call, and its visualization layer, which
  ctdR never used, put Pandoc, cairo, fontconfig, freetype2, libuv and
  glpk among the system requirements of every installation. The
  dependency closure drops from 175 to 116 packages and the system
  requirements from 14 to 7.

* **`maxGSSize` no longer defaults to 500. There is no upper limit.**
  The lower threshold was corrected in this same release because it had
  been inherited from a tool tuned for KEGG and GO; the upper one had
  been inherited from exactly the same place and kept without question.

  The criterion is the one used for the lower threshold, read from the
  other end. A set of *M* genes cannot, even when every one of the *m*
  input genes falls inside it, produce a p-value below
  `choose(M, m) / choose(N, m)`, roughly `(M/N)^m`. That floor rises
  with *M*, so a large enough set is untestable. Whether CTD holds any
  such set is a measurement, not an opinion: with N = 28,571 and an
  input list of 169, the largest still-testable set is about 26,600
  genes, 93% of the universe, while the largest chemical in CTD has
  16,536. No chemical is untestable from above.

  A cap at 500 excluded 265 chemicals, among them benzo(a)pyrene,
  valproic acid, sodium arsenite, bisphenol A, aflatoxin B1 and
  particulate matter: the canonical compounds of toxicology, whose sets
  are large because the literature on them is. This is where CTD parts
  company with GO, in which a large term is one that has stopped meaning
  anything. Removing the cap costs 3% more tests, 8,235 against 7,970.

  In the bundled RNA-seq example the cap excluded dexamethasone, the
  treatment the experiment applied, which without it ranks third of
  8,193 chemicals at an adjusted p-value of 0.0002.

  `maxGSSize` remains settable for anyone who wants a cap of their own.

* The default `minGSSize` for ORA is now **2**, chosen for CTD instead
  of inherited from a general-purpose tool. The previous value, 10,
  came from `enricher()`'s own default, which suits KEGG and GO (median
  set sizes 72 and 11) but not CTD, where the median chemical has 4
  target genes: it silently excluded 7,831 of 11,067 chemicals from
  testing altogether. Coverage goes from 26.8% to 72.0% of chemicals.
  One-gene sets remain excluded on purpose: their hypergeometric
  p-value equals the ratio of input genes to background whichever gene
  they contain, so they measure membership rather than enrichment.

* **Bug fix.** Running an example, knitting the vignette or running
  `R CMD check` overwrote whatever CTD data the user had imported,
  replacing it with the ten-chemical sample those examples run on. The
  examples called `import_CTD()` against the real user cache under
  `tools::R_user_dir()`; only the test suite knew to redirect it. They
  now write to a temporary cache, as does the vignette. Writing outside
  the session temporary directory is also against both CRAN and
  Bioconductor policy, so this was a defect on two counts.

  The examples set `options(ctdR.cache = tempfile())` visibly rather
  than hiding it, because the same option is how a user isolates one
  analysis from their main cache. The test suite now redirects the cache
  from a single `setup.R` as well: five of its six importing files were
  writing to the real one, so running the tests carried the same cost.
  Tests that need a cache of their own restore the suite's, where they
  previously cleared the option, which sent every later test in the run
  back to the user's cache.

* New `ctd_provenance()` returns the record of which CTD release an
  analysis ran on: the `Report created` date CTD stamps into its own file
  header, where the file came from, when it was imported, how many
  chemicals and chemical-gene pairs were retained, and the ctdR version.
  CTD re-releases continuously and does not version its download
  filenames, so that header date is the only thing identifying which
  snapshot a result came from. `import_CTD()` now reads it, caches it and
  prints it, and `enrichment_CTD()` attaches it to every result.

  Where it is stored follows the container: `metadata()` on objects that
  have the slot, which covers the `SummarizedExperiment` GSVA returns and
  the `DataFrame` from importing a `CTDFile`; an attribute on the data
  frames from ORA, GSEA and CAMERA. Use the accessor rather than either
  directly. The record survives subsetting, ordering, `head()` and the
  common dplyr verbs; it does not survive `merge()` or `subset()`, which
  drop attributes, and the accessor says so rather than returning an
  empty answer.

* `import_CTD()` no longer assumes the CTD header is 27 lines long. A CTD
  download has no header row: the field names sit inside the commented
  preamble. The names are now located by finding the commented line that
  lists at least three known CTD field names, taking the last such line,
  and the file is read with `comment = "#"` as suggested in review. The
  hard-coded `skip`, the row dropped afterwards to compensate, and the
  patch that stripped `"# "` from the first column name are all gone.

  This matters beyond tidiness. A fixed line count fails silently: insert
  one comment line upstream and every column shifts, with the analysis
  proceeding on misaligned data. Matching field names fails loudly, and
  it adds no assumption the package was not already making, since those
  names are referenced throughout.

* The bundled sample file now mirrors the structure of a real CTD
  download, header included. It previously carried an uncommented,
  duplicated header row, shaped so the old hard-coded skip would work,
  which meant tests and examples never exercised the format users
  actually have. Its preamble is deliberately a different length from a
  real download's, so that nothing can come to depend on the count again.

* **Bug fix.** The ORA `universe` argument was unusable in the default
  identifier mode, and failed silently. With `gene_id_type = "symbol"`,
  the gene sets are keyed by HGNC symbol and the input gene list is
  converted for that reason, but the universe was passed through as
  given. A universe of Entrez IDs therefore intersected the background
  at nothing, the size filter then removed every gene set, and the call
  returned zero rows instead of an error. The vignette's own example of
  restricting the background to expressed genes shipped in that state.
  The universe is now converted alongside the input; identifiers that
  are not Entrez IDs, and Entrez IDs that do not map, are kept as they
  are, so a universe of symbols or a mixture of the two works too.

* The ORA `universe` argument now accepts any vector of gene
  identifiers, not only a character one. A DE table read back with
  `read.delim()` gives integer Entrez IDs, so the most ordinary use of
  the argument, `universe = de$EntrezID` to restrict the background to
  measured genes, used to fail while the same column passed as the input
  gene list worked. Both are now coerced. The previous backend accepted
  a non-character universe and then silently ignored it, computing
  against a background the caller had not asked for.

* ORA now reports what its size filter removed. A chemical excluded for
  having too few or too many target genes is absent from the results,
  not present with an unremarkable p-value, and the two cases used to be
  indistinguishable. `ora()` emits a message giving how many chemicals
  went untested and on which side of the thresholds they fell.

* **Breaking change.** ORA results now have 10 columns instead of 13.
  `ChemicalID`, `ChemicalName`, `Method`, `PValue`, `PValueAdjusted`,
  `GeneRatio`, `BackgroundRatio`, `EnrichedGenes`, `Count` and
  `FoldEnrichment` are unchanged. Three columns are gone, none of which
  carried information the remaining ones do not:

  - `QValue` held Storey's q-value from the `qvalue` package. On
    result sets of the size a CTD analysis produces it was identical to
    `PValueAdjusted`, because the q-value estimator falls back to
    Benjamini-Hochberg when it cannot estimate the proportion of true
    nulls, so the two columns held the same numbers.
    Reproducing it would mean taking the dependency back for a
    duplicate. Use `PValueAdjusted` for false-discovery control.
  - `RichFactor` and `zScore` came from the previous backend and were
    passed through undocumented: neither appeared in the output schema
    described in the vignette.

  This also fixes a defect. The previous output carried **two** columns
  named `FoldEnrichment`, one from the backend and one from ctdR's own
  rename. Every ORA call raised a duplicated-column warning from
  `merge()`, and `results$FoldEnrichment` returned whichever of the two
  came first. There is now one.

## Internal

* The internal `ora()` engine takes `universe`, `minGSSize` and
  `maxGSSize` as explicit arguments rather than forwarding an opaque
  `...` to another package, so an unrecognised argument now raises an
  error instead of being silently discarded.

* `.parse_ratio()` has been removed. Fold enrichment is computed from
  the counts directly instead of being parsed back out of the
  `"n/d"` strings.

* New `tools/check_internal_params.R`, wired into CI, fails the build
  when a documented function has an argument without its `@param`.
  `R CMD check` skips that cross-check for topics marked
  `\keyword{internal}`, so such an argument used to ship undocumented
  with the check still reporting Status OK. Six internal topics that
  were already in that state have been documented.

* The ORA test suite checks `ora()` against `stats::phyper()` and
  closed-form values rather than against another implementation. After
  the migration the hypergeometric distribution is the reference; an
  equivalence test against `clusterProfiler` would have pinned ctdR's
  correctness to a package it no longer depends on.

# Changes in version 0.99.8

## New features

* `enrichment_CTD()` now accepts a
  `SummarizedExperiment` for the matrix-based methods (`"CAMERA"` and
  `"GSVA"`), in addition to a plain expression matrix. The assay and the
  sample annotation stay in one object, so subsetting or reordering
  samples cannot silently desynchronise them from the group labels used
  to build the design matrix.

* `"GSVA"` returns the container it was given: a matrix in returns a
  numeric matrix, a `SummarizedExperiment` in returns a
  `SummarizedExperiment` whose assay holds the scores and whose
  `colData` is carried over from the input. Per-sample scores therefore
  arrive with the sample annotation needed to interpret them.

* New `assay` argument on `enrichment_CTD()` selects which assay to use
  when the input carries more than one, by name or by index. Defaults to
  the first.

* `plot_CTD()` accepts GSVA scores wrapped in a `SummarizedExperiment`.

## Changes

* The bundled `inst/extdata/GSE311566_subset.rds` example is now stored
  as a `SummarizedExperiment` instead of a `list(expr, coldata)`. Code
  reading it must use `assay(se)` and `se$group` in place of `$expr` and
  `$coldata$group`. `metadata()` records the GEO source, the subsetting
  applied, and the assay units.

* The vignette and the RNA-seq workflow tutorial build a
  `SummarizedExperiment` and carry it through the analysis.

* `SummarizedExperiment` was added to `Imports`. It was already a hard
  transitive dependency through `GSVA`, so the installation footprint is
  unchanged.

# Changes in version 0.99.7

## New features

* New `CTDFile` class (a `BiocIO::BiocFile` subclass) with an `import()`
  method on the `BiocIO` generic. `import(CTDFile(path_or_url))` reads a
  CTD chemical-gene interactions file into a `S4Vectors::DataFrame` of
  validated human interactions, following the Bioconductor import/export
  convention.

* New `as_genesets_CTD()` exports the CTD chemical gene sets as a named
  list keyed by `ChemicalID`, ready for third-party enrichment engines
  such as `EnrichmentBrowser::sbea()`. Supports `id_type` and
  `interaction_types`.

* New `ctd_cache()` retrieves the processed cached tables (`"chemicals"`,
  `"interactions"`) without loading the cache files by hand.

* `import_CTD()` (and `import(CTDFile)`) now accept a local path or a
  URL. Remote URLs are downloaded and cached via `BiocFileCache`; no
  default URL is assumed, and a one-time CTD data-licensing reminder is
  shown on the first remote fetch.

## Changes

* The processed-data cache moved from a `rappdirs` directory of `.rda`
  files to a `BiocFileCache` store under `tools::R_user_dir("ctdR",
  "cache")`. `rappdirs` is no longer a dependency.

* `pAdjustMethod` now accepts any value in `stats::p.adjust.methods`
  (previously restricted to `BH`, `bonferroni`, `fdr`, `none`); its
  documentation references `stats::p.adjust`.

## Documentation

* The vignette gains an "Interoperability with existing Bioconductor
  infrastructure" section (reading via `CTDFile`/`import()`, and handing
  gene sets to `EnrichmentBrowser::sbea()`) plus a note disambiguating
  the "CTD" acronym.

## Internal

* Test coverage raised (package total about 95%); tests exercise the
  `interaction_types` filtering branch and the cache and URL error paths.

# Changes in version 0.99.6

## New features

* `enrichment_CTD()` gains an `interaction_types` argument: a character
  vector of CTD `InteractionActions` values (e.g.
  `"increases^expression"`, `"decreases^expression"`) that filters each
  chemical's gene set at enrichment time. Gene sets are rebuilt on the
  fly from the new `ctd_interactions.rda` cache without re-importing;
  `NULL` (default) retains full backward compatibility. Supported by all
  four methods (ORA, GSEA, CAMERA, GSVA). Requires `import_CTD()` to be
  re-run once to generate `ctd_interactions.rda`.

* `enrichment_CTD()` gains a `gene_id_type` argument (`"symbol"` or
  `"entrez"`) controlling whether the `EnrichedGenes` output column
  reports HGNC symbols (with Entrez ID fallback for unmapped genes) or
  raw Entrez IDs. Default `"symbol"` is backward compatible.

* `enrichment_CTD(method = "ORA")` now forwards `universe`, `minGSSize`,
  and `maxGSSize` to `clusterProfiler::enricher()` via `...`. Setting
  `universe = de$EntrezID` restricts the background to measured genes,
  avoiding inflated fold-enrichment estimates. GSEA and GSVA already
  exposed `minSize`/`maxSize`; CAMERA `...` was already forwarded.

* `import_CTD()` now caches `ctd_interactions.rda`, a long-format table
  of `(ChemicalID, EntrezID, InteractionActions)` triples used by the
  new `interaction_types` filter.

* `import_CTD()` reports elapsed time, chemical count, and unique gene
  count on completion.

* `import_CTD()` detects and warns when the same `ChemicalID` appears
  with multiple `ChemicalName` values (CTD data quality issue; first
  name retained). Reports an informational message when the same name
  is shared by multiple `ChemicalID`s (legitimate parent/derivative
  pairs).

## Bug fixes

* Fixed `EnrichedGenes` column in GSEA output: was incorrectly set to
  data frame row indices instead of gene identifiers. Gene labels are
  now mapped via `AnnotationDbi::mapIds()` in `.run_gsea()`.

* `clusterProfiler::enricher()` messages ("No gene can be mapped",
  "Expected input gene ID", etc.) are now suppressed via
  `suppressMessages()`.

## Documentation

* New pkgdown article `vignettes/articles/tutorial_rnaseq_workflow.Rmd`:
  a complete RNA-seq → chemical enrichment workflow using the full
  GSE311566 dataset (downloaded from GEO at runtime), `limma` DE, all
  four methods with recommended parameters, direction-aware GSEA, and a
  per-method Dexamethasone ranking recap.

* Vignette gains a "Gene set size filters and background universe"
  section documenting `universe`, `minGSSize`/`maxGSSize` (ORA),
  `minSize`/`maxSize` (GSEA, GSVA), and CAMERA's implicit minimum of 2
  genes.

* Added a new `interaction_types` parameter description to the vignette
  explaining the CTD `verb^noun` vocabulary and direction-aware analysis.

* Added a full end-to-end pipeline example at
  `inst/scripts/example_gse311566_full_pipeline.R`. The script
  downloads the complete GSE311566 Female PBMCs normalised-count
  matrix (Dex vs DMSO), computes an a-priori power analysis with
  declared alpha thresholds, runs `limma`-based differential
  expression, and exercises all four enrichment methods (ORA, GSEA,
  CAMERA, GSVA) with BH-adjusted significance cutoffs. It refuses to
  fall back to the bundled toy CTD sample; the user must populate the
  CTD cache with the real chemical-gene interactions file first.
  Complements the vignette (which uses the bundled subset *without*
  alpha cutoffs) by providing the production-shaped example linked
  from the README and the companion paper.

# Changes in version 0.99.5

## Documentation

* Added an end-to-end real-data example to the vignette using a small
  subset of GEO series GSE311566 (human PBMCs, dexamethasone vs.
  vehicle, female donors). The example walks through loading the
  bundled subset, a deliberately minimal base-R differential
  expression with `t.test` + `p.adjust`, and the four ctdR methods
  (ORA, GSEA, CAMERA, GSVA) on the resulting DE.
* Bundled `inst/extdata/GSE311566_subset.rds` (~34 KB) containing
  log2-normalised counts for 1,500 top-variance genes plus the 17
  genes referenced by the toy CTD sample, across 7 samples
  (4 DMSO + 3 Dex).
* Added a reproducible provenance script at
  `inst/scripts/make_gse311566_subset.R` and a per-file documentation
  README at `inst/extdata/README.md`.
* The `enrichment_CTD()` `@examples` block no longer relies on
  `\donttest{}`: CAMERA and GSVA examples now run directly on the
  bundled subset, satisfying the BiocCheck recommendation against
  `\dontrun{}` / `\donttest{}` in man pages.

## Testing

* New `tests/testthat/test-e2e-gse311566.R` runs the full
  data-to-enrichment pipeline on the bundled GSE311566 subset and
  asserts that Dexamethasone (`D003907`) ranks in the top 3 by GSEA
  p-value and in the top 6 by CAMERA p-value, plus structural
  checks on the GSVA output. Guards against silent regressions in
  ID mapping, output schema, or sort order that the vignette and
  man-page examples would only catch as "still runs".

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
