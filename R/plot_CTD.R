#' @title Plot Chemical Enrichment Results
#'
#' @description
#' Visualizes the output of \code{\link{enrichment_CTD}} as a bar plot or dot
#' plot of the top enriched chemicals. The function automatically detects
#' whether the input comes from ORA or GSEA based on column names.
#'
#' @param results A data frame returned by \code{\link{enrichment_CTD}}.
#' @param type Character. Plot type: \code{"bar"} (default) for a bar plot of
#'   fold enrichment colored by adjusted p-value, or \code{"dot"} for a dot
#'   plot with dot size proportional to gene count and color mapped to adjusted
#'   p-value.
#' @param n Integer. Number of top chemicals to display (default 20). Chemicals
#'   are ranked by adjusted p-value.
#' @param title Character. Plot title. If \code{NULL} (default), a title is
#'   generated automatically based on the enrichment method.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' # Import sample data:
#' sample_file <- system.file(
#'     "extdata", "CTD_chem_gene_ixns_sample.csv",
#'     package = "ctdR"
#' )
#' import_CTD(sample_file)
#'
#' genes <- data.frame(
#'     entrez_ids = c("7124", "3569", "7157", "672", "1956"),
#'     pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
#' )
#'
#' ora_results <- enrichment_CTD(genes, method = "ORA")
#' plot_CTD(ora_results, type = "bar")
#' plot_CTD(ora_results, type = "dot")
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_point coord_flip scale_color_gradient
#'   scale_fill_gradient scale_size_continuous labs theme_minimal theme element_text
#' @importFrom stats reorder
#' @export
plot_CTD <- function(results, type = "bar", n = 20, title = NULL) {
    if (!is.data.frame(results) || nrow(results) == 0) {
        stop("'results' must be a non-empty data frame from enrichment_CTD().",
            call. = FALSE)
    }

    type <- match.arg(type, c("bar", "dot"))

    # Detect method: ORA has "Count" and "GeneRatio"; GSEA has "NES" and "ES"
    is_ora <- "Count" %in% colnames(results)

    # Determine padj column
    padj_col <- if ("padj" %in% colnames(results)) "padj" else "pval"

    # Determine gene count column
    if (is_ora) {
        count_col <- "Count"
    } else {
        count_col <- "size"
    }

    # Sort by padj and take top n
    results <- results[order(results[[padj_col]]), ]
    results <- utils::head(results, n)

    # Build plotting data frame
    plot_df <- data.frame(
        ChemicalName = results$ChemicalName,
        foldEnrichment = results$foldEnrichment,
        padj = results[[padj_col]],
        Count = results[[count_col]],
        stringsAsFactors = FALSE
    )

    method_label <- if (is_ora) "ORA" else "GSEA"
    if (is.null(title)) {
        title <- paste("Top", nrow(plot_df),
            "Enriched Chemicals (", method_label, ")")
    }

    if (type == "bar") {
        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(
                x = reorder(ChemicalName, foldEnrichment),
                y = foldEnrichment,
                fill = padj
            )
        ) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::coord_flip() +
            ggplot2::scale_fill_gradient(low = "#D73027", high = "#4575B4",
                name = "Adjusted\np-value") +
            ggplot2::labs(
                x = NULL, y = "Fold Enrichment", title = title
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 10),
                plot.title = ggplot2::element_text(hjust = 0.5)
            )
    } else {
        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(
                x = foldEnrichment,
                y = reorder(ChemicalName, foldEnrichment),
                color = padj,
                size = Count
            )
        ) +
            ggplot2::geom_point() +
            ggplot2::scale_color_gradient(low = "#D73027", high = "#4575B4",
                name = "Adjusted\np-value") +
            ggplot2::scale_size_continuous(
                name = if (is_ora) "Gene Count" else "Gene Set Size"
            ) +
            ggplot2::labs(
                x = "Fold Enrichment", y = NULL, title = title
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 10),
                plot.title = ggplot2::element_text(hjust = 0.5)
            )
    }

    p
}
