#' @title Competitive gene-set enrichment with CAMERA
#'
#' @description
#' Internal function that runs the CAMERA competitive gene-set test
#' (\code{\link[limma]{camera}}) of an expression matrix against the cached
#' CTD chemical gene sets. Unlike ORA/GSEA --- which take a single
#' (ranked) gene list --- CAMERA tests, for each chemical, whether its
#' target genes show a stronger differential signal than the rest of the
#' transcriptome under a user-supplied design and contrast, while accounting
#' for inter-gene correlation.
#'
#' @param expr Numeric expression matrix (genes x samples). \code{rownames(expr)}
#'   must be Entrez IDs or HGNC symbols matching the cached CTD gene sets.
#' @param design Design matrix produced e.g. by \code{\link[stats]{model.matrix}}.
#' @param contrast Either a column number or column name of \code{design}, or
#'   a numeric contrast vector.
#' @param id_type Either \code{"entrez"}, \code{"symbol"}, or \code{NULL} for
#'   auto-detection from \code{rownames(expr)}.
#' @param pAdjustMethod Method passed to \code{\link[stats]{p.adjust}} to
#'   compute the \code{padj} column.
#' @param chemicals_meta Data frame with columns \code{ChemicalID} and
#'   \code{ChemicalName}, used to annotate results.
#' @param cache_dir Directory holding the cached CTD \code{.rda} files.
#' @param ... Forwarded to \code{\link[limma]{camera}} (e.g. \code{inter.gene.cor},
#'   \code{use.ranks}, \code{allow.neg.cor}, \code{trend.var}).
#'
#' @return A data frame with columns \code{ChemicalID}, \code{ChemicalName},
#'   \code{NGenes}, \code{Direction}, \code{Correlation} (if reported),
#'   \code{pvalue}, \code{padj}, sorted by \code{padj} ascending.
#'
#' @keywords internal
.run_camera <- function(expr, design, contrast,
    id_type = NULL,
    pAdjustMethod = "BH",
    chemicals_meta,
    cache_dir,
    ...) {
    .validate_expr_matrix(expr)

    if (is.null(id_type)) {
        id_type <- .detect_id_type(rownames(expr))
    } else {
        id_type <- match.arg(id_type, c("entrez", "symbol"))
    }

    gene_sets <- .load_geneset_list(id_type, cache_dir)

    rn <- as.character(rownames(expr))
    index_list <- lapply(gene_sets, function(g) {
        idx <- which(rn %in% g)
        if (length(idx) < 2L) NULL else idx
    })
    index_list <- index_list[!vapply(index_list, is.null, logical(1))]
    if (length(index_list) == 0L) {
        stop(
            "No chemicals have at least 2 genes matching rownames(x). ",
            "Check that rownames are ", id_type, " IDs.",
            call. = FALSE
        )
    }

    res <- limma::camera(
        y = expr,
        index = index_list,
        design = design,
        contrast = contrast,
        ...
    )

    res$ChemicalID <- rownames(res)
    res$padj <- stats::p.adjust(res$PValue, method = pAdjustMethod)
    res <- merge(res, chemicals_meta, by = "ChemicalID", all.x = TRUE)
    res <- res[order(res$padj), , drop = FALSE]

    colnames(res)[colnames(res) == "PValue"] <- "pvalue"
    keep_fdr <- !colnames(res) %in% "FDR"
    res <- res[, keep_fdr, drop = FALSE]

    front <- c(
        "ChemicalID", "ChemicalName", "NGenes",
        "Direction", "Correlation", "pvalue", "padj"
    )
    front <- intersect(front, colnames(res))
    other <- setdiff(colnames(res), front)
    res <- res[, c(front, other), drop = FALSE]
    rownames(res) <- NULL
    res
}
