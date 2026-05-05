test_that("ora returns enriched results with proper columns", {
    skip_on_cran()
    skip_if_not_installed("clusterProfiler")

    # Build a universe large enough for enricher to find significant results.
    # We need: many terms with many genes, and a gene list that strongly overlaps
    # one term but draws from the broader universe.
    all_genes <- paste0("GENE", 1:500)
    target_genes <- paste0("TARGET", 1:20)

    term2gene <- rbind(
        data.frame(term = "CHEM1", gene = target_genes),
        data.frame(term = "CHEM2", gene = all_genes[1:50]),
        data.frame(term = "CHEM3", gene = all_genes[51:100]),
        data.frame(term = "CHEM4", gene = all_genes[101:200]),
        data.frame(term = "CHEM5", gene = all_genes[201:300])
    )

    # Gene list: all 20 target genes (full overlap with CHEM1) + 10 background
    gene_list <- c(target_genes, all_genes[301:310])

    result <- ora(term2gene, gene_list)

    expect_true(is.data.frame(result))
    expect_true(nrow(result) > 0)
    expect_true("padj" %in% colnames(result))
    expect_true("foldEnrichment" %in% colnames(result))
    expect_true("ID" %in% colnames(result))
    expect_true(is.numeric(result$foldEnrichment))
    expect_true(all(result$foldEnrichment > 0))
})

test_that("ora returns empty data frame when no enrichment found", {
    skip_on_cran()
    skip_if_not_installed("clusterProfiler")

    term2gene <- data.frame(
        term = rep("CHEM1", 2),
        gene = c("TP53", "TNF")
    )
    gene_list <- c("NONEXISTENT_GENE_1", "NONEXISTENT_GENE_2")

    result <- ora(term2gene, gene_list)

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 0)
    expect_true("padj" %in% colnames(result))
    expect_true("foldEnrichment" %in% colnames(result))
})

test_that("ora p.adjust is renamed to padj", {
    skip_on_cran()
    skip_if_not_installed("clusterProfiler")

    all_genes <- paste0("GENE", 1:500)
    target_genes <- paste0("TARGET", 1:20)
    term2gene <- rbind(
        data.frame(term = "CHEM1", gene = target_genes),
        data.frame(term = "CHEM2", gene = all_genes[1:50]),
        data.frame(term = "CHEM3", gene = all_genes[51:100]),
        data.frame(term = "CHEM4", gene = all_genes[101:200])
    )
    gene_list <- c(target_genes, all_genes[201:210])

    result <- ora(term2gene, gene_list)

    # p.adjust should have been renamed
    expect_false("p.adjust" %in% colnames(result))
    expect_true("padj" %in% colnames(result))
})
