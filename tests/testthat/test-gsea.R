.setup_sample_cache_gsea <- function() {
    sample_file <- system.file(
        "extdata", "CTD_chem_gene_ixns_sample.csv",
        package = "ctdR"
    )
    if (!nzchar(sample_file) || !file.exists(sample_file))
        sample_file <- "../../inst/extdata/CTD_chem_gene_ixns_sample.csv"
    if (!file.exists(sample_file)) skip("Sample CTD file not available")
    suppressMessages(suppressWarnings(import_CTD(sample_file)))
    invisible(NULL)
}

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
    # ctdR-added decorations. There is exactly one: the engine used to
    # add a foldEnrichment too, computed as abs(ES) / mean(ES), which
    # described the run rather than the chemical.
    expect_true("Enriched_GENE" %in% colnames(result))
    expect_false("foldEnrichment" %in% colnames(result))
    # ChemicalName / PValueAdjusted / Method are added by the runner;
    # not the engine's job.
    expect_false("ChemicalName" %in% colnames(result))
    expect_false("PValueAdjusted" %in% colnames(result))
    expect_false("Method" %in% colnames(result))
})

test_that("GSEA reports NES and not a fabricated fold enrichment", {
    # gsea() used to add foldEnrichment = abs(ES) / mean(ES). That is not
    # a fold enrichment: the divisor is the mean score across whichever
    # chemicals happened to be tested in the same run, so the value is a
    # property of the run rather than of the chemical, and it duplicates
    # NES, which fgsea computes properly and which is what the field
    # expects. The shared output schema is five columns; method-specific
    # extras differ by method, so GSEA is under no obligation to carry a
    # column ORA has.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    .setup_sample_cache_gsea()
    genes <- data.frame(
        EntrezID = c("7124", "3569", "7157", "672", "1956"),
        pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
    )
    res <- suppressMessages(suppressWarnings(
        enrichment_CTD(genes, method = "GSEA")))

    expect_false("FoldEnrichment" %in% colnames(res))
    expect_true("NormalizedEnrichmentScore" %in% colnames(res))
    # The five shared columns are untouched by the removal.
    expect_true(all(c("ChemicalID", "ChemicalName", "Method",
        "PValue", "PValueAdjusted") %in% colnames(res)))
})
