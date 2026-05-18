.setup_sample_cache_utils <- function() {
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

test_that(".detect_id_type returns 'entrez' for numeric ids", {
    expect_equal(ctdR:::.detect_id_type(c("7124", "3569", "672")),
        "entrez")
})

test_that(".detect_id_type returns 'symbol' for non-numeric ids", {
    expect_equal(ctdR:::.detect_id_type(c("TNF", "IL6", "TP53")),
        "symbol")
})

test_that(".detect_id_type stops on empty / all-NA input", {
    expect_error(
        ctdR:::.detect_id_type(c("", NA_character_, NA_character_)),
        "rownames"
    )
    expect_error(
        ctdR:::.detect_id_type(character(0)),
        "rownames"
    )
})

test_that(".validate_expr_matrix rejects non-matrix input", {
    expect_error(
        ctdR:::.validate_expr_matrix(data.frame(a = 1:3)),
        "numeric matrix"
    )
})

test_that(".validate_expr_matrix rejects matrix without rownames", {
    m <- matrix(rnorm(12), nrow = 3, ncol = 4)
    expect_error(
        ctdR:::.validate_expr_matrix(m),
        "rownames"
    )
})

test_that(".validate_expr_matrix rejects duplicated rownames", {
    m <- matrix(rnorm(12), nrow = 3, ncol = 4,
        dimnames = list(c("A", "A", "B"), paste0("S", 1:4)))
    expect_error(
        ctdR:::.validate_expr_matrix(m),
        "duplicated"
    )
})

test_that(".load_geneset_list returns named list for entrez", {
    skip_on_cran()
    .setup_sample_cache_utils()
    cache_dir <- rappdirs::user_cache_dir("ctdR")
    gs <- ctdR:::.load_geneset_list("entrez", cache_dir)
    expect_true(is.list(gs))
    expect_true(length(gs) > 0)
    expect_true(all(vapply(gs, is.character, logical(1))))
})

test_that(".load_geneset_list returns named list for symbol", {
    skip_on_cran()
    .setup_sample_cache_utils()
    cache_dir <- rappdirs::user_cache_dir("ctdR")
    gs <- ctdR:::.load_geneset_list("symbol", cache_dir)
    expect_true(is.list(gs))
    expect_true(length(gs) > 0)
    expect_true(all(vapply(gs, is.character, logical(1))))
})

test_that("enrichment_CTD works with SYMBOL rownames via id_type='symbol'", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache_utils()
    cache_dir <- rappdirs::user_cache_dir("ctdR")
    gs <- ctdR:::.load_geneset_list("symbol", cache_dir)
    syms <- as.character(unique(unlist(gs)))
    syms <- syms[!is.na(syms) & nzchar(syms)]
    skip_if(length(syms) < 4, "not enough symbol identifiers in sample")

    expr <- matrix(rnorm(length(syms) * 6), nrow = length(syms),
        ncol = 6,
        dimnames = list(syms, paste0("S", 1:6)))
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    res <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2, id_type = "symbol")
    expect_true(is.data.frame(res))
    expect_true(nrow(res) > 0)
})

test_that("enrichment_CTD errors when x is missing", {
    expect_error(enrichment_CTD(), "is required")
})
