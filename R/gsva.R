#' @title Per-sample gene-set scoring with GSVA
#'
#' @description
#' Internal function that runs Gene Set Variation Analysis
#' (\code{\link[GSVA]{gsva}}) on expression data using the cached CTD
#' chemical gene sets. Unlike ORA/GSEA/CAMERA --- which return a single
#' p-value per chemical (group-level inference) --- GSVA produces
#' \strong{per-sample enrichment scores}: one row per chemical, one column
#' per sample. These scores are suitable for downstream clustering,
#' association tests against phenotypes, survival analysis, or heatmap
#' visualization.
#'
#' A \code{SummarizedExperiment} is handed to GSVA unchanged rather than
#' reduced to its assay, because GSVA is itself SE-in/SE-out: the scores
#' then come back in a container that still carries \code{colData}, which
#' is what makes per-sample scores interpretable.
#'
#' @param expr Numeric expression matrix (genes x samples) or a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#'   \code{rownames(expr)} must be Entrez IDs or HGNC symbols matching the
#'   cached CTD gene sets.
#' @param id_type Either \code{"entrez"}, \code{"symbol"}, or \code{NULL} for
#'   auto-detection from \code{rownames(expr)}.
#' @param cache_dir Directory holding the cached CTD \code{.rda} files.
#' @param assay Assay name or index to use when \code{expr} is a
#'   \code{SummarizedExperiment}; \code{NULL} (default) takes the first.
#'   Ignored for a matrix.
#' @param ... Forwarded to \code{\link[GSVA]{gsvaParam}} (e.g. \code{kcdf},
#'   \code{minSize}, \code{maxSize}, \code{tau}, \code{maxDiff}).
#'
#' @return GSVA enrichment scores with CTD chemical IDs in rows and samples
#'   in columns, in the same container as \code{expr}: a numeric matrix for
#'   a matrix input, a \code{SummarizedExperiment} (with \code{colData}
#'   preserved) for a \code{SummarizedExperiment} input.
#'
#' @keywords internal
.run_gsva <- function(expr, id_type = NULL, cache_dir,
    interaction_types = NULL, assay = NULL, ...) {
    ## A SummarizedExperiment is passed through to GSVA rather than reduced
    ## to its assay: GSVA is SE-in/SE-out, so the scores come back in a
    ## container that still carries colData.
    expr <- if (.is_se(expr)) {
        .select_se_assay(expr, assay)
    } else {
        .as_expr_matrix(expr, assay)
    }

    if (is.null(id_type)) {
        id_type <- .detect_id_type(rownames(expr))
    } else {
        id_type <- match.arg(id_type, c("entrez", "symbol"))
    }

    if (!is.null(interaction_types)) {
        gs <- .filter_gene_sets(cache_dir, interaction_types)
        gene_sets <- if (id_type == "entrez") gs$entrez else {
            split(gs$symbols$gene, gs$symbols$term)
        }
    } else {
        gene_sets <- .load_geneset_list(id_type, cache_dir)
    }

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
