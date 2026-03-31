test_that("enrichment_CTD ORA path works with cached data", {
    skip_on_cran()
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("DOSE")
    skip_if_not_installed("AnnotationDbi")
    skip_if_not_installed("org.Hs.eg.db")

    # Set up a temp cache with fake CTD data
    tmp_cache <- file.path(tempdir(), "ctdR_ora_test")
    dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)

    # chemicals lookup
    chemicals <- data.frame(
        ChemicalID = c("CHEM1", "CHEM2"),
        ChemicalName = c("Chemical One", "Chemical Two"),
        stringsAsFactors = FALSE
    )

    # TERM2GENE mapping using real gene symbols
    # Use well-known genes that org.Hs.eg.db can map
    ChemicalName_GeneSymbols <- data.frame(
        term = c(rep("CHEM1", 5), rep("CHEM2", 3)),
        gene = c("TP53", "TNF", "IL6", "BRCA1", "EGFR", "MYC", "KRAS", "AKT1"),
        stringsAsFactors = FALSE
    )

    # Dummy entrez ID list (needed for cache check)
    ChemicalName_GeneEntrezIds <- list(
        CHEM1 = c(7157, 7124, 3569, 672, 1956),
        CHEM2 = c(4609, 3845, 207)
    )

    save(chemicals, file = file.path(tmp_cache, "chemicals.rda"))
    save(ChemicalName_GeneSymbols, file = file.path(tmp_cache, "ChemicalName_GeneSymbols.rda"))
    save(ChemicalName_GeneEntrezIds, file = file.path(tmp_cache, "ChemicalName_GeneEntrezIds.rda"))

    original_cache_fn <- rappdirs::user_cache_dir
    assignInNamespace("user_cache_dir", function(...) tmp_cache, ns = "rappdirs")
    on.exit({
        assignInNamespace("user_cache_dir", original_cache_fn, ns = "rappdirs")
        unlink(tmp_cache, recursive = TRUE)
    })

    # Use Entrez IDs for genes that map to symbols in CHEM1:
    # TP53=7157, TNF=7124, IL6=3569, BRCA1=672, EGFR=1956
    entrez_ids <- data.frame(
        entrez_ids = c("7157", "7124", "3569", "672", "1956"),
        pvalue = c(0.001, 0.002, 0.003, 0.004, 0.005),
        stringsAsFactors = FALSE
    )

    result <- enrichment_CTD(entrez_ids, method = "ORA")

    expect_true(is.data.frame(result))
    # Should have renamed ID -> ChemicalID
    expect_true("ChemicalID" %in% colnames(result))
    # Should have merged ChemicalName
    expect_true("ChemicalName" %in% colnames(result))
    # Description column should be removed
    expect_false("Description" %in% colnames(result))
    # padj and foldEnrichment should exist
    expect_true("padj" %in% colnames(result))
    expect_true("foldEnrichment" %in% colnames(result))
})

test_that("enrichment_CTD GSEA path works with cached data", {
    skip_on_cran()
    skip_if_not_installed("fgsea")

    # Set up a temp cache with fake CTD data
    tmp_cache <- file.path(tempdir(), "ctdR_gsea_test")
    dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)

    chemicals <- data.frame(
        ChemicalID = c("CHEM1", "CHEM2"),
        ChemicalName = c("Chemical One", "Chemical Two"),
        stringsAsFactors = FALSE
    )

    # Pathways as named list of entrez IDs
    ChemicalName_GeneEntrezIds <- list(
        CHEM1 = c(1, 2, 3, 4, 5),
        CHEM2 = c(50, 51, 52, 53, 54)
    )

    save(chemicals, file = file.path(tmp_cache, "chemicals.rda"))
    save(ChemicalName_GeneEntrezIds, file = file.path(tmp_cache, "ChemicalName_GeneEntrezIds.rda"))

    original_cache_fn <- rappdirs::user_cache_dir
    assignInNamespace("user_cache_dir", function(...) tmp_cache, ns = "rappdirs")
    on.exit({
        assignInNamespace("user_cache_dir", original_cache_fn, ns = "rappdirs")
        unlink(tmp_cache, recursive = TRUE)
    })

    set.seed(42)
    entrez_ids <- data.frame(
        entrez_ids = as.character(1:100),
        pvalue = runif(100),
        stringsAsFactors = FALSE
    )

    result <- enrichment_CTD(entrez_ids, method = "GSEA")

    expect_true(is.data.frame(result))
    expect_true("ChemicalID" %in% colnames(result))
    expect_true("ChemicalName" %in% colnames(result))
    expect_true("padj" %in% colnames(result))
    expect_true("foldEnrichment" %in% colnames(result))
    expect_true("Enriched_GENE" %in% colnames(result))
})

test_that("enrichment_CTD GSEA applies pAdjustMethod correctly", {
    skip_on_cran()
    skip_if_not_installed("fgsea")

    tmp_cache <- file.path(tempdir(), "ctdR_gsea_padj_test")
    dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)

    chemicals <- data.frame(
        ChemicalID = c("CHEM1", "CHEM2"),
        ChemicalName = c("Chemical One", "Chemical Two"),
        stringsAsFactors = FALSE
    )

    ChemicalName_GeneEntrezIds <- list(
        CHEM1 = c(1, 2, 3, 4, 5),
        CHEM2 = c(50, 51, 52, 53, 54)
    )

    save(chemicals, file = file.path(tmp_cache, "chemicals.rda"))
    save(ChemicalName_GeneEntrezIds, file = file.path(tmp_cache, "ChemicalName_GeneEntrezIds.rda"))

    original_cache_fn <- rappdirs::user_cache_dir
    assignInNamespace("user_cache_dir", function(...) tmp_cache, ns = "rappdirs")
    on.exit({
        assignInNamespace("user_cache_dir", original_cache_fn, ns = "rappdirs")
        unlink(tmp_cache, recursive = TRUE)
    })

    set.seed(42)
    entrez_ids <- data.frame(
        entrez_ids = as.character(1:100),
        pvalue = runif(100),
        stringsAsFactors = FALSE
    )

    result_bh <- enrichment_CTD(entrez_ids, method = "GSEA", pAdjustMethod = "BH")
    result_bonf <- enrichment_CTD(entrez_ids, method = "GSEA", pAdjustMethod = "bonferroni")

    expect_true(is.data.frame(result_bh))
    expect_true(is.data.frame(result_bonf))
    # Bonferroni is more conservative, so padj values should be >= BH values
    # Merge by ChemicalID to compare
    merged <- merge(result_bh[, c("ChemicalID", "padj")],
        result_bonf[, c("ChemicalID", "padj")],
        by = "ChemicalID", suffixes = c("_bh", "_bonf")
    )
    if (nrow(merged) > 0) {
        expect_true(all(merged$padj_bonf >= merged$padj_bh))
    }
})

test_that("enrichment_CTD GSEA handles NA entrez_ids", {
    skip_on_cran()
    skip_if_not_installed("fgsea")

    tmp_cache <- file.path(tempdir(), "ctdR_gsea_na_test")
    dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)

    chemicals <- data.frame(
        ChemicalID = c("CHEM1"),
        ChemicalName = c("Chemical One"),
        stringsAsFactors = FALSE
    )

    ChemicalName_GeneEntrezIds <- list(
        CHEM1 = c(1, 2, 3, 4, 5)
    )

    save(chemicals, file = file.path(tmp_cache, "chemicals.rda"))
    save(ChemicalName_GeneEntrezIds, file = file.path(tmp_cache, "ChemicalName_GeneEntrezIds.rda"))

    original_cache_fn <- rappdirs::user_cache_dir
    assignInNamespace("user_cache_dir", function(...) tmp_cache, ns = "rappdirs")
    on.exit({
        assignInNamespace("user_cache_dir", original_cache_fn, ns = "rappdirs")
        unlink(tmp_cache, recursive = TRUE)
    })

    set.seed(42)
    entrez_ids <- data.frame(
        entrez_ids = c(as.character(1:50), NA, NA),
        pvalue = c(runif(50), 0.5, 0.5),
        stringsAsFactors = FALSE
    )

    # Should not error — NAs are removed internally
    result <- enrichment_CTD(entrez_ids, method = "GSEA")
    expect_true(is.data.frame(result))
})
