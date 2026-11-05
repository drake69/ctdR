#' @title Gene Set Enrichment Analysis (GSEA)
#'
#' @description
#' Internal engine that performs Gene Set Enrichment Analysis via
#' \code{\link[fgsea]{fgsea}}. Detects chemicals whose known target
#' genes cluster toward the extremes of a ranked gene list.
#'
#' This function returns the raw GSEA result with two GSEA-specific
#' decorations (\code{Enriched_GENE}) and
#' harmonized column names (\code{ChemicalID}, \code{pvalue}). The
#' adjusted-p, chemical-name join, column ordering and sort are
#' applied downstream by
#' \code{\link{.format_enrichment_result}} so that the schema stays
#' consistent across ORA / GSEA / CAMERA.
#'
#' @param ChemicalName_GeneEntrezIds A named list where each element
#'   is a character vector of Entrez gene IDs associated with a CTD
#'   chemical. Names are CTD chemical IDs.
#' @param gene_table A data frame with at least two columns:
#'   \describe{
#'     \item{\code{EntrezID}}{Entrez gene IDs (character).}
#'     \item{(second column)}{P-values, used as the
#'       \code{-log10(value)} fallback ranking when no \code{stat} column
#'       is present. Values outside [0, 1] stop the run: a fold change
#'       there would be ranked upside down, and its negative values
#'       dropped.}
#'     \item{\code{stat} (optional)}{Numeric ranking statistic used
#'       directly by \code{fgsea} (higher = ranked first). Recommended:
#'       the moderated t-statistic from \code{limma::eBayes()} or
#'       \code{-log10(p) * sign(logFC)}. When present, suppresses the
#'       directionality warning. A \code{stat} column whose values all
#'       lie in [0, 1] raises a warning, because p-values put there would
#'       rank the least significant genes first.}
#'   }
#'
#' @param ... Additional arguments forwarded to
#'   \code{\link[fgsea]{fgseaMultilevel}} (e.g. \code{nproc} to set the
#'   number of parallel workers, \code{eps} for p-value precision). With
#'   the \code{-log10(p)} fallback every value is positive, so
#'   \code{scoreType} defaults to \code{"pos"} unless given here.
#'
#' @return A data frame with \code{fgsea::fgseaMultilevel}'s native columns
#'   (\code{pathway}, \code{pval}, \code{ES}, \code{NES}, \code{size},
#'   \code{leadingEdge}) plus two ctdR-added columns
#'   (\code{Enriched_GENE}). Not yet sorted and not
#'   yet joined with chemical metadata; pass to
#'   \code{\link{.format_enrichment_result}} for the canonical shape.
#'
#' @keywords internal
gsea <- function(ChemicalName_GeneEntrezIds, gene_table, ...) {
    fgsea_args <- list(...)
    if ("stat" %in% colnames(gene_table)) {
        stats <- gene_table$stat
        if (.looks_like_pvalues(stats))
            warning(
                "GSEA: every value in the 'stat' column lies in [0, 1], ",
                "as p-values do. fgsea ranks the highest value first, so ",
                "p-values there put the least significant genes at the top. ",
                "Use -log10(p), -log10(p) * sign(logFC) or a t-statistic.",
                call. = FALSE
            )
    } else {
        if (!.looks_like_pvalues(gene_table[, 2]))
            stop(
                "GSEA: without a 'stat' column the second column must hold ",
                "p-values in [0, 1], ranked as -log10(p). Column '",
                colnames(gene_table)[2], "' has values outside [0, 1]: ",
                "a fold change there would be ranked upside down and its ",
                "negative values dropped. Add a 'stat' column with the ",
                "ranking statistic, or put the p-values second.",
                call. = FALSE
            )
        warning(
            "GSEA: no 'stat' column found in input. Genes are ranked by ",
            "-log10(second column), without sign, and scored one-sided: the ",
            "test asks whether a chemical's genes change more than the ",
            "others, in either direction. To ask about one direction, add a ",
            "signed statistic (e.g. the moderated t-statistic from ",
            "limma::eBayes()) as a column named 'stat'.",
            call. = FALSE
        )
        stats <- -log10(gene_table[, 2] + 1e-10)
        # Every -log10(p) is positive: the genes at the bottom of the list
        # are the unchanged ones, not the opposite tail, so only the top
        # can carry an enrichment and the test is one-sided.
        if (is.null(fgsea_args$scoreType)) fgsea_args$scoreType <- "pos"
    }
    names(stats) <- gene_table$EntrezID
    stats <- sort(stats, decreasing = TRUE)

    fgsea_results <- do.call(fgsea::fgseaMultilevel, c(
        list(pathways = ChemicalName_GeneEntrezIds, stats = stats),
        fgsea_args
    ))
    fgsea_results <- as.data.frame(fgsea_results)

    # fgsea calls its primary key "pathway" (its API targets biological
    # pathways), but in ctdR those IDs are CTD chemical IDs. Lift to the
    # semantically correct name here so callers never see "pathway".
    names(fgsea_results)[
        names(fgsea_results) == "pathway"
    ] <- "ChemicalID"

    fgsea_results$Enriched_GENE <- .annotate_genes(
        fgsea_results, gene_table, ChemicalName_GeneEntrezIds
    )
    fgsea_results
}

#' Do these values look like p-values?
#'
#' GSEA ranks the highest value first, so a p-value column used as is
#' turns the ranking upside down. Checking the range is the one test that
#' needs no knowledge of where the column came from.
#' @param values Numeric vector; \code{NA} values are ignored.
#' @return \code{TRUE} when \code{values} is numeric and every non-missing
#'   value lies in [0, 1], \code{FALSE} otherwise.
#' @keywords internal
.looks_like_pvalues <- function(values) {
    observed <- values[!is.na(values)]
    is.numeric(values) && length(observed) > 0 &&
        all(observed >= 0 & observed <= 1)
}

#' Annotate GSEA results with enriched gene symbols
#' @param results Data frame of enrichment results carrying a
#'   \code{ChemicalID} column.
#' @param gene_table Data frame with an \code{EntrezID} column and a
#'   \code{GeneLabel} column holding the label to report.
#' @param chemical_sets Named list mapping ChemicalID to the Entrez IDs
#'   of that chemical's target genes.
#' @return A character vector of comma-separated gene symbols per chemical.
#' @keywords internal
.annotate_genes <- function(results, gene_table, chemical_sets) {
    vapply(
        results$ChemicalID,
        function(x) {
            matched <- gene_table$EntrezID %in% chemical_sets[[x]]
            paste(gene_table[matched, "GeneLabel"], collapse = ", ")
        },
        character(1)
    )
}
