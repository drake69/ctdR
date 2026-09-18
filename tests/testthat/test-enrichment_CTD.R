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

test_that("universe is refused loudly by the methods that cannot use it", {
    # Only ORA takes a background it cannot infer. GSEA ranks the whole
    # list supplied, and CAMERA and GSVA intersect the gene sets with
    # rownames(x), so those three already have the measured genes as
    # their background. Accepting the argument and dropping it would
    # leave a caller believing they had narrowed something.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    genes <- data.frame(EntrezID = c("7124", "3569"), pvalue = c(0.001, 0.003))
    se <- readRDS(system.file("extdata", "GSE311566_subset.rds",
        package = "ctdR"))
    design <- stats::model.matrix(~ se$group)

    expect_warning(
        suppressMessages(enrichment_CTD(genes, method = "GSEA",
            universe = "7124")),
        "ignored for \"GSEA\"")
    expect_warning(
        suppressMessages(enrichment_CTD(se, method = "CAMERA",
            design = design, contrast = 2, universe = "7124")),
        "rownames\\(x\\)")
    expect_warning(
        suppressMessages(enrichment_CTD(se, method = "GSVA",
            universe = "7124")),
        "ignored for \"GSVA\"")

    # ORA accepts it without comment: that is where it belongs.
    expect_no_warning(
        suppressMessages(enrichment_CTD(genes, method = "ORA",
            universe = c("7124", "3569", "7157"))))
})

test_that("alpha derives the gene list and the background from one table", {
    # The failure this prevents: filtering first and passing a separate
    # universe asks the caller to reconnect two things that were together
    # a moment earlier. Given alpha, they cannot come apart.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    de <- data.frame(
        EntrezID = c("7124", "3569", "7157", "672", "1956", "836", "7422"),
        padj = c(0.001, 0.003, 0.01, 0.02, 0.30, 0.40, 0.50)
    )

    from_alpha <- suppressMessages(
        enrichment_CTD(de, method = "ORA", alpha = 0.05))
    by_hand <- suppressMessages(enrichment_CTD(
        de[de$padj < 0.05, ], method = "ORA", universe = de$EntrezID))

    expect_equal(from_alpha, by_hand)
    expect_gt(nrow(from_alpha), 0)
})

test_that("alpha reports the selection and the full table as background", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    de <- data.frame(
        EntrezID = c("7124", "3569", "7157", "672", "1956", "836", "7422"),
        padj = c(0.001, 0.003, 0.01, 0.02, 0.30, 0.40, 0.50)
    )

    # The background is every row, the selected genes included. The
    # hypergeometric draws n genes from an urn of N and the drawn ones
    # were in the urn; describing the background as the complement would
    # leave nothing to test.
    expect_message(enrichment_CTD(de, method = "ORA", alpha = 0.05),
        "4 genes selected")
    expect_message(enrichment_CTD(de, method = "ORA", alpha = 0.05),
        "background = all 7 rows of the table")

    # And the reported N is the whole table, not the 3 non-selected rows.
    res <- suppressMessages(enrichment_CTD(de, method = "ORA", alpha = 0.05))
    expect_identical(sub(".*/", "", res$BackgroundRatio[1]), "7")
})

test_that("alpha and universe are mutually exclusive", {
    # Not a precedence rule applied in silence: with alpha the background
    # is already decided, so universe has nothing left to say and passing
    # both means the caller believes something untrue.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    de <- data.frame(EntrezID = c("7124", "3569"), padj = c(0.001, 0.30))

    expect_error(
        enrichment_CTD(de, method = "ORA", alpha = 0.05, universe = "7124"),
        "either 'alpha' or 'universe'")
})

test_that("alpha_column decides which p-value the threshold judges", {
    # A real result table has several numeric columns and the second is
    # usually a fold change: limma::topTable() puts logFC there.
    # Thresholding a position instead of a name would filter on the wrong
    # quantity without saying so.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    de <- data.frame(
        EntrezID = c("7124", "3569", "7157", "672", "1956", "836", "7422"),
        log2FC = c(2.1, -1.8, 1.5, -2.2, 0.3, 0.1, -0.2),
        pvalue = c(0.0001, 0.001, 0.004, 0.01, 0.2, 0.3, 0.4),
        padj = c(0.001, 0.003, 0.01, 0.02, 0.30, 0.40, 0.50)
    )

    by_name <- suppressMessages(enrichment_CTD(de, method = "ORA",
        alpha = 0.05, alpha_column = "padj"))
    by_index <- suppressMessages(enrichment_CTD(de, method = "ORA",
        alpha = 0.05, alpha_column = 4))
    expect_equal(by_name, by_index)

    # A stricter threshold on a different column selects differently.
    expect_message(enrichment_CTD(de, method = "ORA", alpha = 0.005,
        alpha_column = "pvalue"), "3 genes selected")

    expect_error(enrichment_CTD(de, method = "ORA", alpha = 0.05,
        alpha_column = "fdr"), "not in 'x'")
    expect_error(enrichment_CTD(de, method = "ORA", alpha = 0.05,
        alpha_column = 9), "out of range")
})

test_that("alpha refuses a table it cannot threshold", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    expect_error(
        suppressMessages(enrichment_CTD(
            data.frame(EntrezID = c("7124", "3569"), label = c("a", "b")),
            method = "ORA", alpha = 0.05)),
        "not numeric")
    expect_error(
        suppressMessages(enrichment_CTD(
            data.frame(EntrezID = "7124", padj = 0.9),
            method = "ORA", alpha = 0.05)),
        "No gene is below alpha")
})

test_that("a missing EntrezID column is named, not left to AnnotationDbi", {
    # The README shipped an example using `entrez_ids`, which failed with
    # "mapIds must have at least one key to match against": an error from
    # a package the caller never invoked, naming neither the column nor
    # the function that wanted it.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache()
    wrong <- data.frame(entrez_ids = c("7124", "3569"), pvalue = c(0.001, 0.003))

    for (m in c("ORA", "GSEA")) {
        expect_error(enrichment_CTD(wrong, method = m), "needs a column named")
        expect_error(enrichment_CTD(wrong, method = m), "entrez_ids")
    }
})
