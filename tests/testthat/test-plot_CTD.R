test_that("plot_CTD bar plot works with ORA results", {
    skip_if_not_installed("ggplot2")

    ora_df <- data.frame(
        ChemicalID = paste0("D00", 1:5),
        ChemicalName = paste("Chemical", LETTERS[1:5]),
        GeneRatio = c("3/10", "2/10", "4/10", "1/10", "5/10"),
        BgRatio = c("20/500", "30/500", "15/500", "40/500", "10/500"),
        pvalue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        padj = c(0.005, 0.01, 0.03, 0.05, 0.1),
        foldEnrichment = c(3.5, 2.1, 4.2, 1.3, 5.0),
        geneID = c("A/B/C", "D/E", "F/G/H/I", "J", "K/L/M/N/O"),
        Count = c(3L, 2L, 4L, 1L, 5L),
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(ora_df, type = "bar")

    expect_s3_class(p, "ggplot")
})

test_that("plot_CTD dot plot works with ORA results", {
    skip_if_not_installed("ggplot2")

    ora_df <- data.frame(
        ChemicalID = paste0("D00", 1:5),
        ChemicalName = paste("Chemical", LETTERS[1:5]),
        GeneRatio = c("3/10", "2/10", "4/10", "1/10", "5/10"),
        BgRatio = c("20/500", "30/500", "15/500", "40/500", "10/500"),
        pvalue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        padj = c(0.005, 0.01, 0.03, 0.05, 0.1),
        foldEnrichment = c(3.5, 2.1, 4.2, 1.3, 5.0),
        geneID = c("A/B/C", "D/E", "F/G/H/I", "J", "K/L/M/N/O"),
        Count = c(3L, 2L, 4L, 1L, 5L),
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(ora_df, type = "dot")

    expect_s3_class(p, "ggplot")
})

test_that("plot_CTD works with GSEA results", {
    skip_if_not_installed("ggplot2")

    gsea_df <- data.frame(
        ChemicalID = paste0("D00", 1:3),
        ChemicalName = paste("Chemical", LETTERS[1:3]),
        pval = c(0.001, 0.01, 0.05),
        padj = c(0.005, 0.03, 0.1),
        ES = c(0.65, 0.45, 0.30),
        NES = c(1.8, 1.5, 1.2),
        size = c(15L, 20L, 10L),
        leadingEdge = I(list(c("A", "B"), c("C", "D"), c("E"))),
        foldEnrichment = c(2.5, 1.8, 1.2),
        Enriched_GENE = c("A, B", "C, D", "E"),
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(gsea_df, type = "bar")

    expect_s3_class(p, "ggplot")
})

test_that("plot_CTD limits to n top results", {
    skip_if_not_installed("ggplot2")

    ora_df <- data.frame(
        ChemicalID = paste0("D", 1:30),
        ChemicalName = paste("Chemical", 1:30),
        pvalue = seq(0.001, 0.03, length.out = 30),
        padj = seq(0.005, 0.1, length.out = 30),
        foldEnrichment = seq(5, 1, length.out = 30),
        Count = rep(3L, 30),
        geneID = rep("A/B/C", 30),
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(ora_df, n = 10)

    expect_s3_class(p, "ggplot")
    # The plot data should have only 10 rows
    expect_equal(nrow(p$data), 10)
})

test_that("plot_CTD errors on empty data frame", {
    expect_error(
        plot_CTD(data.frame()),
        "non-empty data frame"
    )
})

test_that("plot_CTD accepts custom title", {
    skip_if_not_installed("ggplot2")

    ora_df <- data.frame(
        ChemicalID = "D001",
        ChemicalName = "Chemical A",
        pvalue = 0.001,
        padj = 0.005,
        foldEnrichment = 3.5,
        Count = 3L,
        geneID = "A/B/C",
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(ora_df, title = "My Custom Title")

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$title, "My Custom Title")
})
