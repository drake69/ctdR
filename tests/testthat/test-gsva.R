.setup_sample_cache_gsva <- function() {
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

.synthetic_expr_gsva <- function(seed = 7, n_samples = 6) {
    cache_dir <- rappdirs::user_cache_dir("ctdR")
    e <- new.env(parent = emptyenv())
    load(file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"), envir = e)
    ids <- as.character(unique(unlist(e$ChemicalName_GeneEntrezIds)))
    set.seed(seed)
    matrix(
        rnorm(length(ids) * n_samples),
        nrow = length(ids), ncol = n_samples,
        dimnames = list(ids, paste0("S", seq_len(n_samples)))
    )
}

test_that("GSVA returns a numeric matrix with chemicals as rows and samples as columns", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    .setup_sample_cache_gsva()

    expr <- .synthetic_expr_gsva()
    suppressMessages({
        scores <- enrichment_CTD(expr, method = "GSVA")
    })

    expect_true(is.matrix(scores))
    expect_true(is.numeric(scores))
    expect_equal(ncol(scores), ncol(expr))
    expect_equal(colnames(scores), colnames(expr))
    expect_true(!is.null(rownames(scores)))
    expect_true(all(nchar(rownames(scores)) > 0))
})

test_that("GSVA score range is reasonable for the GSVA method", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    .setup_sample_cache_gsva()

    expr <- .synthetic_expr_gsva()
    suppressMessages({
        scores <- enrichment_CTD(expr, method = "GSVA")
    })

    expect_true(all(scores >= -1.0001 & scores <= 1.0001))
})

test_that("GSVA errors when x is a data.frame", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    .setup_sample_cache_gsva()

    df <- data.frame(EntrezID = "7124", value = 0.01)
    expect_error(
        enrichment_CTD(df, method = "GSVA"),
        "numeric matrix"
    )
})

test_that("GSVA errors when no gene set has >= 2 matched genes", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    .setup_sample_cache_gsva()

    expr <- matrix(rnorm(20 * 6), nrow = 20, ncol = 6,
        dimnames = list(
            paste0("UNMATCHED_", seq_len(20)),
            paste0("S", 1:6)
        ))

    expect_error(
        suppressMessages(enrichment_CTD(expr, method = "GSVA")),
        "matching rownames"
    )
})

test_that("GSVA respects explicit id_type override", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    .setup_sample_cache_gsva()

    expr <- .synthetic_expr_gsva()
    suppressMessages({
        s_auto <- enrichment_CTD(expr, method = "GSVA")
        s_explicit <- enrichment_CTD(expr, method = "GSVA",
            id_type = "entrez")
    })
    expect_equal(dim(s_auto), dim(s_explicit))
    expect_equal(rownames(s_auto), rownames(s_explicit))
})

test_that("GSVA errors when x has no rownames", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    .setup_sample_cache_gsva()

    expr <- matrix(rnorm(30), nrow = 5, ncol = 6)
    expect_error(
        suppressMessages(enrichment_CTD(expr, method = "GSVA")),
        "rownames"
    )
})

test_that("GSVA passes additional arguments through to gsvaParam", {
    skip_on_cran()
    skip_if_not_installed("GSVA")
    .setup_sample_cache_gsva()

    expr <- .synthetic_expr_gsva()
    suppressMessages({
        s_def <- enrichment_CTD(expr, method = "GSVA")
        # minSize = 5 should filter out most/all chemicals in the small sample
        s_strict <- enrichment_CTD(expr, method = "GSVA", minSize = 5)
    })
    expect_true(nrow(s_strict) <= nrow(s_def))
})
