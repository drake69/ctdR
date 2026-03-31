test_that("enrichment_CTD errors when CTD data not imported", {
    # Temporarily override user_cache_dir to an empty temp dir
    tmp_cache <- file.path(tempdir(), "ctdR_empty_cache")
    dir.create(tmp_cache, showWarnings = FALSE)
    original_cache_fn <- rappdirs::user_cache_dir
    assignInNamespace("user_cache_dir", function(...) tmp_cache, ns = "rappdirs")
    on.exit({
        assignInNamespace("user_cache_dir", original_cache_fn, ns = "rappdirs")
        unlink(tmp_cache, recursive = TRUE)
    })

    expect_error(
        enrichment_CTD(data.frame(entrez_ids = "7124", value = 0.01)),
        "CTD data not found"
    )
    expect_error(
        enrichment_CTD(data.frame(entrez_ids = "7124", value = 0.01)),
        "import_CTD"
    )
    expect_error(
        enrichment_CTD(data.frame(entrez_ids = "7124", value = 0.01)),
        "CTD_chem_gene_ixns.csv.gz"
    )
})

test_that("enrichment_CTD errors on invalid pAdjustMethod", {
    expect_error(
        enrichment_CTD(data.frame(entrez_ids = "7124", value = 0.01),
            pAdjustMethod = "invalid"
        ),
        "pAdjustMethod"
    )
})

test_that("enrichment_CTD error mentions download URL", {
    tmp_cache <- file.path(tempdir(), "ctdR_empty_cache2")
    dir.create(tmp_cache, showWarnings = FALSE)
    original_cache_fn <- rappdirs::user_cache_dir
    assignInNamespace("user_cache_dir", function(...) tmp_cache, ns = "rappdirs")
    on.exit({
        assignInNamespace("user_cache_dir", original_cache_fn, ns = "rappdirs")
        unlink(tmp_cache, recursive = TRUE)
    })

    expect_error(
        enrichment_CTD(data.frame(entrez_ids = "7124", value = 0.01)),
        "ctdbase.org"
    )
})
