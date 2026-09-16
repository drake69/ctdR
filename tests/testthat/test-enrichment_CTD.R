test_that("enrichment_CTD errors when CTD data not imported", {
    # Temporarily override user_cache_dir to an empty temp dir
    tmp_cache <- file.path(tempdir(), "ctdR_empty_cache")
    dir.create(tmp_cache, showWarnings = FALSE)
    options(ctdR.cache = tmp_cache)
    on.exit({
        .restore_ctd_cache()
        unlink(tmp_cache, recursive = TRUE)
    })

    expect_error(
        enrichment_CTD(data.frame(EntrezID = "7124", value = 0.01)),
        "CTD data not found"
    )
    expect_error(
        enrichment_CTD(data.frame(EntrezID = "7124", value = 0.01)),
        "import_CTD"
    )
    expect_error(
        enrichment_CTD(data.frame(EntrezID = "7124", value = 0.01)),
        "CTD_chem_gene_ixns.csv.gz"
    )
})

test_that("enrichment_CTD errors on invalid pAdjustMethod", {
    expect_error(
        enrichment_CTD(data.frame(EntrezID = "7124", value = 0.01),
            pAdjustMethod = "invalid"
        ),
        "pAdjustMethod"
    )
})

test_that("pAdjustMethod accepts every stats::p.adjust.methods value", {
    df <- data.frame(EntrezID = "7124", value = 0.01)
    empty <- file.path(tempdir(), "ctdR_padj_validate")
    dir.create(empty, showWarnings = FALSE)
    on.exit(unlink(empty, recursive = TRUE))

    # A valid p.adjust method passes the pAdjustMethod check and only then
    # fails on the (empty) cache -- proving the method itself was accepted.
    for (m in stats::p.adjust.methods) {
        expect_error(
            ctdR:::.validate_enrichment_args(df, "ORA", NULL, NULL, m, empty),
            "CTD data not found"
        )
    }
    expect_error(
        ctdR:::.validate_enrichment_args(df, "ORA", NULL, NULL, "nope", empty),
        "must be one of"
    )
})

test_that("enrichment_CTD error mentions download URL", {
    tmp_cache <- file.path(tempdir(), "ctdR_empty_cache2")
    dir.create(tmp_cache, showWarnings = FALSE)
    options(ctdR.cache = tmp_cache)
    on.exit({
        .restore_ctd_cache()
        unlink(tmp_cache, recursive = TRUE)
    })

    expect_error(
        enrichment_CTD(data.frame(EntrezID = "7124", value = 0.01)),
        "ctdbase.org"
    )
})

.setup_sample_cache <- function() {
    sample_file <- system.file(
        "extdata", "CTD_chem_gene_ixns_sample.csv",
        package = "ctdR"
    )
    if (!nzchar(sample_file) || !file.exists(sample_file)) {
        sample_file <- "../../inst/extdata/CTD_chem_gene_ixns_sample.csv"
    }
    if (!file.exists(sample_file)) {
        skip("Sample CTD file not available")
    }
    suppressMessages(suppressWarnings(import_CTD(sample_file)))
    invisible(NULL)
}

test_that("an Entrez universe works in the default symbol mode", {
    # The gene sets are keyed by symbol by default and the input list is
    # converted for that reason; the universe has to make the same trip.
    # Left as Entrez IDs it intersects the background at nothing, and an
    # empty background yields an empty result rather than an error, which
    # is how the vignette shipped an example returning zero rows.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    genes <- data.frame(
        EntrezID = c("7124", "3569", "7157", "672", "1956"),
        pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
    )

    entrez_universe <- suppressMessages(enrichment_CTD(genes, method = "ORA",
        universe = c(genes$EntrezID, "7422", "836")))
    expect_gt(nrow(entrez_universe), 0)

    # A universe already given as symbols must keep working: mapping a
    # vector with no valid Entrez key errors in AnnotationDbi.
    symbol_universe <- suppressMessages(enrichment_CTD(genes, method = "ORA",
        universe = c("TNF", "IL6", "TP53", "CASP3", "AKT1")))
    expect_gt(nrow(symbol_universe), 0)

    # And a mixture of the two, which is what a hand-assembled list is.
    mixed <- suppressMessages(enrichment_CTD(genes, method = "ORA",
        universe = c("7124", "IL6", "7157")))
    expect_gt(nrow(mixed), 0)

    # Narrowing the universe must actually narrow the background.
    wide <- suppressMessages(enrichment_CTD(genes, method = "ORA"))
    expect_lte(
        as.integer(sub(".*/", "", symbol_universe$BackgroundRatio[1])),
        as.integer(sub(".*/", "", wide$BackgroundRatio[1]))
    )
})

test_that("an Entrez universe still works in entrez mode", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    genes <- data.frame(
        EntrezID = c("7124", "3569", "7157"),
        pvalue = c(0.001, 0.003, 0.01)
    )
    res <- suppressMessages(enrichment_CTD(genes, method = "ORA",
        gene_id_type = "entrez", universe = c(genes$EntrezID, "7422")))

    expect_gt(nrow(res), 0)
})
