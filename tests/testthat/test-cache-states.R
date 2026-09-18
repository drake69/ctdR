# The two states a user is ever in: no CTD data imported, or some. The
# example script branches on exactly this, and until now it decided by
# rebuilding a path that the package had stopped using, so it could be
# wrong in both directions at once. These tests pin the contract the
# script relies on: what the package says when the cache is empty, and
# what it hands back when it is not.
#
# The bundled sample is a ten-chemical toy, which is the point: the
# states are what is being tested, not the biology.

.empty_cache <- function() {
    dir <- file.path(tempdir(), paste0("ctdR_empty_", basename(tempfile())))
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    dir
}

.sample_csv <- function() {
    f <- system.file("extdata", "CTD_chem_gene_ixns_sample.csv",
        package = "ctdR")
    if (!nzchar(f) || !file.exists(f))
        f <- "../../inst/extdata/CTD_chem_gene_ixns_sample.csv"
    if (!file.exists(f)) skip("Sample CTD file not available")
    f
}

test_that("with an empty cache, every entry point says so and stops", {
    skip_on_cran()

    dir <- .empty_cache()
    options(ctdR.cache = dir)
    on.exit({
        .restore_ctd_cache()
        unlink(dir, recursive = TRUE)
    }, add = TRUE)

    # The analysis refuses, and the message says what to do about it.
    expect_error(
        enrichment_CTD(data.frame(EntrezID = "7124", pvalue = 0.01),
            method = "ORA"),
        "CTD data not found")
    expect_error(
        enrichment_CTD(data.frame(EntrezID = "7124", pvalue = 0.01),
            method = "ORA"),
        "import_CTD")

    # The cache accessor names the resource it could not find and where
    # it looked. The example script turns this error into its own.
    expect_error(ctd_cache("chemicals"), "not found")
    expect_error(ctd_cache("chemicals"), "import_CTD")

    # Provenance degrades to a warning and NULL, never to a silent empty
    # answer: a blank line in a methods section is worse than being told.
    expect_warning(prov <- ctd_provenance(), "No CTD provenance in the cache")
    expect_null(prov)

    # And the gene sets cannot be exported either.
    expect_error(as_genesets_CTD("entrez"), "not found")
})

test_that("with a populated cache, the script's two questions are answerable", {
    # Step E of the example script asks the package for the cached
    # chemicals and for the release date. If either return shape changes,
    # the script breaks silently, so both are pinned here.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    dir <- .empty_cache()
    options(ctdR.cache = dir)
    on.exit({
        .restore_ctd_cache()
        unlink(dir, recursive = TRUE)
    }, add = TRUE)

    suppressMessages(suppressWarnings(import_CTD(.sample_csv())))

    chemicals <- ctd_cache("chemicals")
    expect_s3_class(chemicals, "data.frame")
    expect_true(all(c("ChemicalID", "ChemicalName") %in% colnames(chemicals)))
    expect_gt(nrow(chemicals), 0)

    prov <- ctd_provenance()
    expect_s3_class(prov, "ctd_provenance")
    expect_false(is.na(prov$report_created))
    expect_identical(prov$n_chemicals, nrow(chemicals))

    # And the analysis runs.
    res <- suppressMessages(enrichment_CTD(
        data.frame(EntrezID = c("7124", "3569"), pvalue = c(0.001, 0.003)),
        method = "ORA", universe = c("7124", "3569", "7157", "672")))
    expect_s3_class(res, "data.frame")
    expect_gt(nrow(res), 0)
})

test_that("the bundled sample is small enough for the script to refuse it", {
    # The example script stops when the cache holds fewer than 1,000
    # chemicals, so that it cannot report "0 significant" from a toy
    # universe and have it read as a finding. That guard is only
    # meaningful if the toy really is that small.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    dir <- .empty_cache()
    options(ctdR.cache = dir)
    on.exit({
        .restore_ctd_cache()
        unlink(dir, recursive = TRUE)
    }, add = TRUE)

    suppressMessages(suppressWarnings(import_CTD(.sample_csv())))
    expect_lt(nrow(ctd_cache("chemicals")), 1000L)
})

test_that("importing twice replaces the cache rather than accumulating", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    dir <- .empty_cache()
    options(ctdR.cache = dir)
    on.exit({
        .restore_ctd_cache()
        unlink(dir, recursive = TRUE)
    }, add = TRUE)

    suppressMessages(suppressWarnings(import_CTD(.sample_csv())))
    first <- ctd_provenance()
    n_first <- nrow(ctd_cache("chemicals"))

    suppressMessages(suppressWarnings(import_CTD(.sample_csv())))
    second <- ctd_provenance()

    expect_identical(nrow(ctd_cache("chemicals")), n_first)
    expect_identical(second$report_created, first$report_created)
})
