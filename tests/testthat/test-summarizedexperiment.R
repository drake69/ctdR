## SummarizedExperiment input for the matrix-based methods (CAMERA, GSVA).
##
## The contract under test: a SummarizedExperiment and its assay must be
## interchangeable as input and produce identical numbers, and GSVA must
## return the container it was given, with colData intact.

.setup_sample_cache_se <- function() {
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
    suppressMessages(import_CTD(sample_file))
    invisible(NULL)
}

.synthetic_expr_se <- function(seed = 42, n_samples = 6) {
    ids <- as.character(unique(unlist(as_genesets_CTD("entrez"))))
    set.seed(seed)
    matrix(
        rnorm(length(ids) * n_samples),
        nrow = length(ids), ncol = n_samples,
        dimnames = list(ids, paste0("S", seq_len(n_samples)))
    )
}

.make_se <- function(expr, extra_assay = FALSE) {
    skip_if_not_installed("SummarizedExperiment")
    assays <- list(logcounts = expr)
    if (extra_assay) {
        ## Decoy in first position: any test that selects "logcounts" by
        ## name would still pass if selection were ignored, unless the
        ## default (first assay) differs from it.
        assays <- list(raw = expr * 10, logcounts = expr)
    }
    SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        colData = S4Vectors::DataFrame(
            group = factor(rep(c("ctrl", "treat"), each = ncol(expr) / 2)),
            row.names = colnames(expr)
        )
    )
}

## --- helpers ---------------------------------------------------------------

test_that(".is_se recognises a SummarizedExperiment and nothing else", {
    skip_if_not_installed("SummarizedExperiment")
    expr <- matrix(1, nrow = 2, ncol = 2,
        dimnames = list(c("a", "b"), c("s1", "s2")))
    expect_true(ctdR:::.is_se(.make_se(expr)))
    expect_false(ctdR:::.is_se(expr))
    expect_false(ctdR:::.is_se(data.frame(a = 1)))
    expect_false(ctdR:::.is_se(NULL))
})

test_that(".assay_index resolves NULL, index and name", {
    a <- list(first = 1, second = 2)
    expect_equal(ctdR:::.assay_index(a, NULL), 1L)
    expect_equal(ctdR:::.assay_index(a, 2), 2L)
    expect_equal(ctdR:::.assay_index(a, "second"), 2L)
})

test_that(".assay_index rejects bad selectors with a usable message", {
    a <- list(first = 1, second = 2)
    expect_error(ctdR:::.assay_index(a, 3), "out of range")
    expect_error(ctdR:::.assay_index(a, 0), "out of range")
    expect_error(ctdR:::.assay_index(a, "missing"), "not found")
    expect_error(ctdR:::.assay_index(a, c(1, 2)), "single assay")
    expect_error(ctdR:::.assay_index(list(), 1), "no assays")
})

test_that(".as_expr_matrix extracts the requested assay", {
    expr <- .synthetic_expr_se()
    se <- .make_se(expr, extra_assay = TRUE)

    ## default is the first assay, here the decoy
    expect_equal(ctdR:::.as_expr_matrix(se), expr * 10)
    expect_equal(ctdR:::.as_expr_matrix(se, "logcounts"), expr)
    expect_equal(ctdR:::.as_expr_matrix(se, 2), expr)
    ## a plain matrix passes through untouched
    expect_identical(ctdR:::.as_expr_matrix(expr), expr)
})

test_that(".as_expr_matrix inherits the matrix rowname validation", {
    skip_if_not_installed("SummarizedExperiment")
    bad <- matrix(1, nrow = 2, ncol = 2,
        dimnames = list(c("a", "a"), c("s1", "s2")))
    expect_error(ctdR:::.as_expr_matrix(.make_se(bad)), "duplicated rownames")
})

test_that(".select_se_assay keeps the container and one assay", {
    expr <- .synthetic_expr_se()
    se <- .make_se(expr, extra_assay = TRUE)
    out <- ctdR:::.select_se_assay(se, "logcounts")

    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(length(SummarizedExperiment::assays(out)), 1L)
    expect_equal(as.matrix(SummarizedExperiment::assay(out)), expr)
    expect_identical(out$group, se$group)
})

## --- CAMERA ----------------------------------------------------------------

test_that("CAMERA gives identical results for a matrix and its SE", {
    skip_on_cran()
    skip_if_not_installed("limma")
    skip_if_not_installed("SummarizedExperiment")
    .setup_sample_cache_se()

    expr <- .synthetic_expr_se()
    se <- .make_se(expr)
    d <- model.matrix(~ se$group)

    from_matrix <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2)
    from_se <- enrichment_CTD(se, method = "CAMERA",
        design = d, contrast = 2)

    expect_s3_class(from_se, "data.frame")
    expect_equal(from_se, from_matrix)
})

test_that("CAMERA honours the assay selector", {
    skip_on_cran()
    skip_if_not_installed("limma")
    skip_if_not_installed("SummarizedExperiment")
    .setup_sample_cache_se()

    expr <- .synthetic_expr_se()
    se <- .make_se(expr, extra_assay = TRUE)
    d <- model.matrix(~ se$group)

    picked <- enrichment_CTD(se, method = "CAMERA", design = d,
        contrast = 2, assay = "logcounts")
    expected <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2)
    expect_equal(picked, expected)
})

## --- GSVA ------------------------------------------------------------------

test_that("GSVA returns the container it was given", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    skip_if_not_installed("SummarizedExperiment")
    .setup_sample_cache_se()

    expr <- .synthetic_expr_se()
    se <- .make_se(expr)

    suppressMessages({
        from_matrix <- enrichment_CTD(expr, method = "GSVA")
        from_se <- enrichment_CTD(se, method = "GSVA")
    })

    expect_true(is.matrix(from_matrix))
    expect_s4_class(from_se, "SummarizedExperiment")
    ## Same scores either way. Two attributes ride on the matrix output and
    ## not on the assay of an SE: "geneSets", which GSVA adds, and the CTD
    ## provenance, which goes into metadata() when there is a slot for it.
    ## Neither is a score, so the comparison strips both.
    scores <- from_matrix
    attr(scores, "geneSets") <- NULL
    attr(scores, "ctd_provenance") <- NULL
    expect_equal(
        as.matrix(SummarizedExperiment::assay(from_se)),
        scores
    )
    ## and the provenance is on both, in the place each container provides
    expect_s3_class(ctd_provenance(from_matrix), "ctd_provenance")
    expect_s3_class(ctd_provenance(from_se), "ctd_provenance")
    expect_identical(ctd_provenance(from_matrix)$report_created,
        ctd_provenance(from_se)$report_created)

    ## and the sample annotation is still attached
    expect_identical(from_se$group, se$group)
    expect_identical(colnames(from_se), colnames(se))
})

## --- validation ------------------------------------------------------------

test_that("enrichment_CTD rejects a data frame for the matrix methods", {
    skip_on_cran()
    .setup_sample_cache_se()
    expect_error(
        enrichment_CTD(data.frame(a = 1:3), method = "GSVA"),
        "numeric matrix .* or a SummarizedExperiment"
    )
})

## --- plotting --------------------------------------------------------------

test_that("plot_CTD accepts GSVA scores wrapped in a SummarizedExperiment", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    skip_if_not_installed("SummarizedExperiment")
    .setup_sample_cache_se()

    expr <- .synthetic_expr_se()
    se <- .make_se(expr)
    suppressMessages({
        scores_se <- enrichment_CTD(se, method = "GSVA")
        scores_mat <- enrichment_CTD(expr, method = "GSVA")
    })

    p_se <- plot_CTD(scores_se)
    p_mat <- plot_CTD(scores_mat)
    expect_s3_class(p_se, "ggplot")
    expect_equal(p_se$data, p_mat$data)
})

test_that("plot_CTD error message names every accepted input", {
    expect_error(plot_CTD(data.frame()), "SummarizedExperiment")
})
