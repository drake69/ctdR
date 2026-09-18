.setup_sample_cache_plot <- function() {
    f <- system.file("extdata", "CTD_chem_gene_ixns_sample.csv",
        package = "ctdR")
    if (!nzchar(f) || !file.exists(f))
        f <- "../../inst/extdata/CTD_chem_gene_ixns_sample.csv"
    if (!file.exists(f)) skip("Sample CTD file not available")
    suppressMessages(suppressWarnings(import_CTD(f)))
    invisible(NULL)
}

test_that("plot_CTD bar plot works with ORA results", {
    skip_if_not_installed("ggplot2")

    ora_df <- data.frame(
        ChemicalID = paste0("D00", 1:5),
        ChemicalName = paste("Chemical", LETTERS[1:5]),
        GeneRatio = c("3/10", "2/10", "4/10", "1/10", "5/10"),
        BackgroundRatio = c("20/500", "30/500", "15/500", "40/500", "10/500"),
        PValue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        PValueAdjusted = c(0.005, 0.01, 0.03, 0.05, 0.1),
        FoldEnrichment = c(3.5, 2.1, 4.2, 1.3, 5.0),
        EnrichedGenes = c("A/B/C", "D/E", "F/G/H/I", "J", "K/L/M/N/O"),
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
        BackgroundRatio = c("20/500", "30/500", "15/500", "40/500", "10/500"),
        PValue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        PValueAdjusted = c(0.005, 0.01, 0.03, 0.05, 0.1),
        FoldEnrichment = c(3.5, 2.1, 4.2, 1.3, 5.0),
        EnrichedGenes = c("A/B/C", "D/E", "F/G/H/I", "J", "K/L/M/N/O"),
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
        PValue = c(0.001, 0.01, 0.05),
        PValueAdjusted = c(0.005, 0.03, 0.1),
        EnrichmentScore = c(0.65, 0.45, 0.30),
        NormalizedEnrichmentScore = c(1.8, 1.5, 1.2),
        GeneSetSize = c(15L, 20L, 10L),
        LeadingEdge = I(list(c("A", "B"), c("C", "D"), c("E"))),
        FoldEnrichment = c(2.5, 1.8, 1.2),
        EnrichedGenes = c("A, B", "C, D", "E"),
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
        PValue = seq(0.001, 0.03, length.out = 30),
        PValueAdjusted = seq(0.005, 0.1, length.out = 30),
        FoldEnrichment = seq(5, 1, length.out = 30),
        Count = rep(3L, 30),
        EnrichedGenes = rep("A/B/C", 30),
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
        GeneSetSize = c(15L, 20L, 10L, 8L, 12L),
        Direction = c("Up", "Down", "Up", "Down", "Up"),
        PValue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        PValueAdjusted = c(0.005, 0.01, 0.03, 0.05, 0.1),
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
        GeneSetSize = c(15L, 20L, 10L, 8L, 12L),
        Direction = c("Up", "Down", "Up", "Down", "Up"),
        PValue = c(0.001, 0.005, 0.01, 0.02, 0.05),
        PValueAdjusted = c(0.005, 0.01, 0.03, 0.05, 0.1),
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
        PValue = 0.001,
        PValueAdjusted = 0.005,
        FoldEnrichment = 3.5,
        Count = 3L,
        EnrichedGenes = "A/B/C",
        stringsAsFactors = FALSE
    )

    p <- plot_CTD(ora_df, title = "My Custom Title")

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$title, "My Custom Title")
})

test_that("each method is plotted on its own measure of effect", {
    # plot_CTD() drew FoldEnrichment for both ORA and GSEA. When the
    # fabricated GSEA fold was removed, the function broke at the point
    # where the pipeline draws its figures: a contract change whose
    # consumers had not been checked. Each method now plots the quantity
    # it actually has, and the axis says which.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("ggplot2")

    .setup_sample_cache_plot()
    genes <- data.frame(
        EntrezID = c("7124", "3569", "7157", "672", "1956"),
        pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
    )

    ora <- suppressMessages(enrichment_CTD(genes, method = "ORA"))
    gsea <- suppressMessages(suppressWarnings(
        enrichment_CTD(genes, method = "GSEA")))

    expect_identical(plot_CTD(ora, type = "bar")$labels$y, "Fold Enrichment")
    expect_identical(plot_CTD(ora, type = "dot")$labels$x, "Fold Enrichment")
    expect_identical(plot_CTD(gsea, type = "bar")$labels$y,
        "Normalized Enrichment Score (NES)")
    expect_identical(plot_CTD(gsea, type = "dot")$labels$x,
        "Normalized Enrichment Score (NES)")

    # And a frame missing the column says which one, rather than failing
    # inside data.frame() with a row-count mismatch.
    expect_error(plot_CTD(ora[, setdiff(colnames(ora), "FoldEnrichment")]),
        "FoldEnrichment' is missing")
})
