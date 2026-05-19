#' @title Per-sample gene-set scoring with GSVA
#'
#' @description
#' Internal function that runs Gene Set Variation Analysis
#' (\code{\link[GSVA]{gsva}}) on an expression matrix using the cached CTD
#' chemical gene sets. Unlike ORA/GSEA/CAMERA --- which return a single
#' p-value per chemical (group-level inference) --- GSVA produces a matrix
#' of \strong{per-sample enrichment scores}: one row per chemical, one column
#' per sample. These scores are suitable for downstream clustering,
#' association tests against phenotypes, survival analysis, or heatmap
#' visualization.
#'
#' @param expr Numeric expression matrix (genes x samples).
#'   \code{rownames(expr)} must be Entrez IDs or HGNC symbols matching the
#'   cached CTD gene sets.
#' @param id_type Either \code{"entrez"}, \code{"symbol"}, or \code{NULL} for
#'   auto-detection from \code{rownames(expr)}.
#' @param cache_dir Directory holding the cached CTD \code{.rda} files.
#' @param ... Forwarded to \code{\link[GSVA]{gsvaParam}} (e.g. \code{kcdf},
#'   \code{minSize}, \code{maxSize}, \code{tau}, \code{maxDiff}).
#'
#' @return A numeric matrix of GSVA enrichment scores with CTD chemical IDs
#'   in rows and samples in columns.
#'
#' @keywords internal
.run_gsva <- function(expr, id_type = NULL, cache_dir, ...) {
    .validate_expr_matrix(expr)

    if (is.null(id_type)) {
        id_type <- .detect_id_type(rownames(expr))
    } else {
        id_type <- match.arg(id_type, c("entrez", "symbol"))
    }

    gene_sets <- .load_geneset_list(id_type, cache_dir)

    rn <- as.character(rownames(expr))
    gene_sets <- lapply(gene_sets, function(g) intersect(as.character(g), rn))
    gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= 2L]
    if (length(gene_sets) == 0L) {
        stop(
            "No chemicals have at least 2 genes matching rownames(x). ",
            "Check that rownames are ", id_type, " IDs.",
            call. = FALSE
        )
    }

    param <- GSVA::gsvaParam(
        exprData = expr,
        geneSets = gene_sets,
        ...
    )

    GSVA::gsva(param)
}
