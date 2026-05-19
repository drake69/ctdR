test_that("gsea engine returns fgsea-native + ctdR decorations", {
    skip_on_cran()
    skip_if_not_installed("fgsea")

    chemical_sets <- list(
        CHEM1 = c(1, 2, 3),
        CHEM2 = c(4, 5, 6)
    )

    set.seed(42)
    gene_table <- data.frame(
        EntrezID = as.character(1:100),
        pvalue = runif(100)
    )

    result <- gsea(chemical_sets, gene_table)

    expect_true(is.data.frame(result))
    # gsea() lifts fgsea's "pathway" column to the semantically
    # correct "ChemicalID" and keeps fgsea's native casing for the
    # rest (pval, ES, NES, size, leadingEdge). The PascalCase rename
    # happens downstream in .format_enrichment_result.
    expect_true("ChemicalID" %in% colnames(result))
    expect_false("pathway" %in% colnames(result))
    expect_true("pval" %in% colnames(result))
    expect_true("ES" %in% colnames(result))
    expect_true("NES" %in% colnames(result))
    expect_true("size" %in% colnames(result))
    expect_true("leadingEdge" %in% colnames(result))
    # ctdR-added decorations
    expect_true("foldEnrichment" %in% colnames(result))
    expect_true("Enriched_GENE" %in% colnames(result))
    # ChemicalName / PValueAdjusted / Method are added by the runner;
    # not the engine's job.
    expect_false("ChemicalName" %in% colnames(result))
    expect_false("PValueAdjusted" %in% colnames(result))
    expect_false("Method" %in% colnames(result))
})
