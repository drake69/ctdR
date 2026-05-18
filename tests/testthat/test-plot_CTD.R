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
        "non-empty"
    )
})

test_that("plot_CTD bar plot works with CAMERA results", {
    skip_if_not_installed("ggplot2")

    camera_df <- data.frame(
        ChemicalID = paste0("D00", 1:5),
        ChemicalName = paste("Chemical", LETTERS[1:5]),
        NGenes = c(15L, 20L, 10L, 8L, 12L),
        Direction = c("Up", "Down", "Up", "Down", "Up"),
        pvalue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        padj = c(0.005, 0.01, 0.03, 0.05, 0.1),
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(camera_df, type = "bar")
    expect_s3_class(p, "ggplot")
    # CAMERA branch uses negLog10padj on x axis
    expect_true("negLog10padj" %in% colnames(p$data))
})

test_that("plot_CTD dot plot works with CAMERA results", {
    skip_if_not_installed("ggplot2")

    camera_df <- data.frame(
        ChemicalID = paste0("D00", 1:5),
        ChemicalName = paste("Chemical", LETTERS[1:5]),
        NGenes = c(15L, 20L, 10L, 8L, 12L),
        Direction = c("Up", "Down", "Up", "Down", "Up"),
        pvalue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        padj = c(0.005, 0.01, 0.03, 0.05, 0.1),
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(camera_df, type = "dot")
    expect_s3_class(p, "ggplot")
})

test_that("plot_CTD heatmap works with GSVA matrix", {
    skip_if_not_installed("ggplot2")

    set.seed(1)
    scores <- matrix(rnorm(30, sd = 0.4), nrow = 5, ncol = 6,
        dimnames = list(paste0("D00", 1:5), paste0("S", 1:6)))

    p <- plot_CTD(scores)
    expect_s3_class(p, "ggplot")
    # GSVA branch: data has Score column
    expect_true("Score" %in% colnames(p$data))
})

test_that("plot_CTD limits GSVA heatmap to top-n by variance", {
    skip_if_not_installed("ggplot2")

    set.seed(2)
    scores <- matrix(rnorm(80), nrow = 20, ncol = 4,
        dimnames = list(paste0("D", 1:20), paste0("S", 1:4)))
    # boost the variance of chemicals 1 and 2
    scores[1, ] <- c(-1, -1, 1, 1)
    scores[2, ] <- c(1, -1, -1, 1)

    p <- plot_CTD(scores, n = 3)
    # 3 chemicals x 4 samples = 12 tiles
    expect_equal(nrow(p$data), 3L * 4L)
})

test_that("plot_CTD errors on empty numeric matrix", {
    expect_error(
        plot_CTD(matrix(numeric(0), nrow = 0, ncol = 0)),
        "rows|empty|non-empty"
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
