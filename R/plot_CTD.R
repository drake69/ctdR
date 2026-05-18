#' @title Plot Chemical Enrichment Results
#'
#' @description
#' Visualizes the output of \code{\link{enrichment_CTD}}. The function
#' auto-detects the enrichment method from the input class and column
#' structure:
#' \itemize{
#'   \item \strong{ORA} / \strong{GSEA} (data frame with \code{Count} or
#'     \code{NES}): bar or dot plot of fold enrichment.
#'   \item \strong{CAMERA} (data frame with \code{Direction} and \code{NGenes}):
#'     bar or dot plot of \eqn{-\log_{10}(\mathrm{padj})}, colored by
#'     direction of enrichment.
#'   \item \strong{GSVA} (numeric matrix \code{chemicals x samples}): heatmap
#'     of per-sample enrichment scores for the top-N chemicals selected by
#'     score variance across samples.
#' }
#'
#' @param results A data frame or numeric matrix returned by
#'   \code{\link{enrichment_CTD}}.
#' @param type Character. Plot type for tabular results:
#'   \code{"bar"} (default) or \code{"dot"}. Ignored for GSVA matrix input.
#' @param n Integer. Number of top chemicals to display (default 20).
#'   Selection criterion: ascending \code{padj} for tabular results, descending
#'   score variance across samples for GSVA matrices.
#' @param title Character. Plot title. If \code{NULL} (default), a title is
#'   generated automatically based on the detected method.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' sample_file <- system.file(
#'     "extdata", "CTD_chem_gene_ixns_sample.csv",
#'     package = "ctdR"
#' )
#' import_CTD(sample_file)
#'
#' # ORA / GSEA
#' genes <- data.frame(
#'     entrez_ids = c("7124", "3569", "7157", "672", "1956"),
#'     pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
#' )
#' ora_results <- enrichment_CTD(genes, method = "ORA")
#' plot_CTD(ora_results, type = "bar")
#' plot_CTD(ora_results, type = "dot")
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_tile coord_flip
#' @importFrom ggplot2 scale_color_gradient scale_color_gradient2
#' @importFrom ggplot2 scale_fill_gradient scale_fill_gradient2
#' @importFrom ggplot2 scale_size_continuous scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text element_blank
#' @importFrom stats reorder var
#' @export
plot_CTD <- function(results, type = "bar", n = 20, title = NULL) {
    if (is.matrix(results) && is.numeric(results)) {
        return(.plot_gsva_heatmap(results, n = n, title = title))
    }

    if (!is.data.frame(results) || nrow(results) == 0) {
        stop("'results' must be a non-empty data frame or numeric matrix ",
            "from enrichment_CTD().",
            call. = FALSE)
    }

    type <- match.arg(type, c("bar", "dot"))

    is_camera <- all(c("Direction", "NGenes") %in% colnames(results))
    is_ora <- "Count" %in% colnames(results)

    if (is_camera) {
        return(.plot_camera(results, type = type, n = n, title = title))
    }

    padj_col <- if ("padj" %in% colnames(results)) "padj" else "pval"
    count_col <- if (is_ora) "Count" else "size"

    results <- results[order(results[[padj_col]]), ]
    results <- utils::head(results, n)

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

#' Plot CAMERA enrichment results
#' @keywords internal
.plot_camera <- function(results, type, n, title) {
    results <- results[order(results$padj), ]
    results <- utils::head(results, n)

    plot_df <- data.frame(
        ChemicalName = results$ChemicalName,
        negLog10padj = -log10(pmax(results$padj, .Machine$double.eps)),
        NGenes = results$NGenes,
        Direction = factor(results$Direction, levels = c("Up", "Down")),
        stringsAsFactors = FALSE
    )

    if (is.null(title)) {
        title <- paste("Top", nrow(plot_df),
            "Enriched Chemicals ( CAMERA )")
    }

    dir_colors <- c(Up = "#D73027", Down = "#4575B4")

    if (type == "bar") {
        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(
                x = reorder(ChemicalName, negLog10padj),
                y = negLog10padj,
                fill = Direction
            )
        ) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::coord_flip() +
            ggplot2::scale_fill_manual(values = dir_colors,
                name = "Direction") +
            ggplot2::labs(
                x = NULL,
                y = expression(-log[10] * "(adjusted p-value)"),
                title = title
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
                x = negLog10padj,
                y = reorder(ChemicalName, negLog10padj),
                color = Direction,
                size = NGenes
            )
        ) +
            ggplot2::geom_point() +
            ggplot2::scale_color_manual(values = dir_colors,
                name = "Direction") +
            ggplot2::scale_size_continuous(name = "Gene Set Size") +
            ggplot2::labs(
                x = expression(-log[10] * "(adjusted p-value)"),
                y = NULL,
                title = title
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 10),
                plot.title = ggplot2::element_text(hjust = 0.5)
            )
    }
    p
}

#' Plot GSVA scores as a sample x chemical heatmap
#' @keywords internal
.plot_gsva_heatmap <- function(scores, n, title) {
    if (nrow(scores) == 0L || ncol(scores) == 0L) {
        stop("'results' is an empty matrix; nothing to plot.",
            call. = FALSE)
    }
    vars <- apply(scores, 1, stats::var)
    keep <- order(vars, decreasing = TRUE)[seq_len(min(n, nrow(scores)))]
    mat <- scores[keep, , drop = FALSE]

    long <- data.frame(
        Chemical = factor(rep(rownames(mat), ncol(mat)),
            levels = rev(rownames(mat))),
        Sample = factor(rep(colnames(mat), each = nrow(mat)),
            levels = colnames(mat)),
        Score = as.numeric(mat),
        stringsAsFactors = FALSE
    )

    if (is.null(title)) {
        title <- paste("GSVA Scores (top",
            nrow(mat), "chemicals by variance )")
    }

    p <- ggplot2::ggplot(
        long,
        ggplot2::aes(x = Sample, y = Chemical, fill = Score)
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
            low = "#4575B4", mid = "white", high = "#D73027",
            midpoint = 0, name = "GSVA\nscore"
        ) +
        ggplot2::labs(x = NULL, y = NULL, title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            axis.text.y = ggplot2::element_text(size = 9),
            plot.title = ggplot2::element_text(hjust = 0.5),
            panel.grid = ggplot2::element_blank()
        )
    p
}
