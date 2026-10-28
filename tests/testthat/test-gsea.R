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

test_that("GSEA refuses a second column that is not a p-value", {
    # Without 'stat' the second column is ranked as -log10(value). A fold
    # change there, as limma::topTable() puts it, would be ranked upside
    # down and its negative values would become NaN and vanish in sort().
    skip_on_cran()
    skip_if_not_installed("fgsea")

    chemical_sets <- list(CHEM1 = c("1", "2", "3"))
    gene_table <- data.frame(
        EntrezID = as.character(1:20),
        logFC = seq(-2, 2, length.out = 20)
    )
    expect_error(gsea(chemical_sets, gene_table), "outside \\[0, 1\\]")
})

test_that("GSEA warns when the stat column holds p-values", {
    # fgsea ranks the highest value first: p-values used as 'stat' put
    # the least significant genes at the top.
    skip_on_cran()
    skip_if_not_installed("fgsea")

    set.seed(42)
    chemical_sets <- list(CHEM1 = c("1", "2", "3"))
    gene_table <- data.frame(
        EntrezID = as.character(1:50),
        pvalue = runif(50),
        stat = runif(50)
    )
    # fgsea adds its own warning on an all-positive 'stat'.
    expect_warning(
        expect_warning(gsea(chemical_sets, gene_table), "lies in \\[0, 1\\]"),
        "scoreType"
    )
})

test_that("GSEA scores the -log10(p) fallback one-sided", {
    # Every -log10(p) is positive, so fgsea's default two-sided score
    # does not apply: the fallback passes scoreType = "pos", and a value
    # given by the caller still wins.
    skip_on_cran()
    skip_if_not_installed("fgsea")

    set.seed(42)
    chemical_sets <- list(CHEM1 = as.character(1:5),
        CHEM2 = as.character(30:40))
    gene_table <- data.frame(
        EntrezID = as.character(1:100),
        pvalue = sort(runif(100))
    )
    fallback_warnings <- character(0)
    withCallingHandlers(
        result <- gsea(chemical_sets, gene_table),
        warning = function(condition) {
            fallback_warnings <<- c(fallback_warnings,
                conditionMessage(condition))
            invokeRestart("muffleWarning")
        }
    )
    expect_false(any(grepl("scoreType", fallback_warnings)))
    expect_true(all(result$ES >= 0))

    expect_warning(
        expect_warning(gsea(chemical_sets, gene_table, scoreType = "std"),
            "scoreType"),
        "no 'stat' column"
    )
})
