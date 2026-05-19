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
    suppressMessages(import_CTD(sample_file))
    invisible(NULL)
}

.synthetic_expr <- function(seed = 42, n_samples = 6) {
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

test_that("CAMERA returns a data frame with expected columns", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache()

    expr <- .synthetic_expr()
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    res <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2)

    expect_true(is.data.frame(res))
    expect_true(all(c("ChemicalID", "ChemicalName", "GeneSetSize",
        "Direction", "PValue", "PValueAdjusted", "Method") %in%
        colnames(res)))
    expect_equal(unique(res$Method), "CAMERA")
    expect_true(nrow(res) > 0)
})

test_that("CAMERA results are sorted by PValueAdjusted ascending", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache()

    expr <- .synthetic_expr()
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    res <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2)
    if (nrow(res) > 1) {
        expect_true(all(diff(res$PValueAdjusted) >= 0))
    }
})

test_that("CAMERA PValueAdjusted matches the requested adjustment method", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache()

    expr <- .synthetic_expr()
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    res_bh <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2, pAdjustMethod = "BH")
    res_bf <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2, pAdjustMethod = "bonferroni")
    res_none <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2, pAdjustMethod = "none")

    expect_equal(res_none$PValueAdjusted, res_none$PValue)
    expect_true(all(
        res_bf$PValueAdjusted >= res_bh$PValueAdjusted - 1e-12
    ))
})

test_that("CAMERA errors when design or contrast missing", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache()

    expr <- .synthetic_expr()
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    expect_error(
        enrichment_CTD(expr, method = "CAMERA"),
        "'design' and 'contrast' are required"
    )
    expect_error(
        enrichment_CTD(expr, method = "CAMERA", design = d),
        "'design' and 'contrast' are required"
    )
})

test_that("CAMERA errors when x is a data.frame", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache()

    df <- data.frame(EntrezID = "7124", value = 0.01)
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    expect_error(
        enrichment_CTD(df, method = "CAMERA",
            design = d, contrast = 2),
        "numeric matrix"
    )
})

test_that("ORA errors when x is a matrix", {
    skip_on_cran()
    .setup_sample_cache()

    expr <- .synthetic_expr()
    expect_error(
        enrichment_CTD(expr, method = "ORA"),
        "data.frame"
    )
})

test_that("CAMERA accepts both numeric- and character-typed rownames as Entrez", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache()

    expr <- .synthetic_expr()
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    res_auto <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2)
    res_explicit <- enrichment_CTD(expr, method = "CAMERA",
        design = d, contrast = 2, id_type = "entrez")

    expect_equal(res_auto[, c("ChemicalID", "PValue")],
        res_explicit[, c("ChemicalID", "PValue")])
})

test_that("CAMERA errors when no gene set has >=2 matched genes", {
    skip_on_cran()
    skip_if_not_installed("limma")
    .setup_sample_cache()

    # rownames that don't match any chemical gene set
    expr <- matrix(rnorm(20 * 6), nrow = 20, ncol = 6,
        dimnames = list(
            paste0("UNMATCHED_", seq_len(20)),
            paste0("S", 1:6)
        ))
    grp <- factor(rep(c("ctrl", "treat"), each = 3))
    d <- model.matrix(~ grp)

    expect_error(
        enrichment_CTD(expr, method = "CAMERA",
            design = d, contrast = 2),
        "matching rownames"
    )
})

