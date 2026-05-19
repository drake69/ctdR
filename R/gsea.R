#' @title Gene Set Enrichment Analysis (GSEA)
#'
#' @description
#' Internal engine that performs Gene Set Enrichment Analysis via
#' \code{\link[fgsea]{fgsea}}. Detects chemicals whose known target
#' genes cluster toward the extremes of a ranked gene list.
#'
#' This function returns the raw GSEA result with two GSEA-specific
#' decorations (\code{foldEnrichment}, \code{Enriched_GENE}) and
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
#'     \item{(second column)}{Numeric values used for ranking (e.g.
#'       p-values). Converted internally to \code{-log10(value)} for
#'       the ranking statistic.}
#'   }
#'
#' @return A data frame with \code{fgsea::fgsea}'s native columns
#'   (\code{pathway}, \code{pval}, \code{ES}, \code{NES}, \code{size},
#'   \code{leadingEdge}, \code{nMoreExtreme}) plus two ctdR-added
#'   columns (\code{foldEnrichment}, \code{Enriched_GENE}). Not yet
#'   sorted and not yet joined with chemical metadata; pass to
#'   \code{\link{.format_enrichment_result}} for the canonical shape.
#'
#' @keywords internal
gsea <- function(ChemicalName_GeneEntrezIds, gene_table) {
    stats <- -log10(gene_table[, 2] + 1e-10)
    names(stats) <- gene_table$EntrezID
    stats <- sort(stats, decreasing = TRUE)

    fgsea_results <- fgsea::fgsea(
        pathways = ChemicalName_GeneEntrezIds, stats = stats
    )
    fgsea_results <- as.data.frame(fgsea_results)

    # fgsea calls its primary key "pathway" (its API targets biological
    # pathways), but in ctdR those IDs are CTD chemical IDs. Lift to the
    # semantically correct name here so callers never see "pathway".
    names(fgsea_results)[
        names(fgsea_results) == "pathway"
    ] <- "ChemicalID"

    fgsea_results$foldEnrichment <- (
        abs(fgsea_results$ES) / mean(fgsea_results$ES)
    )
    gene_table$Symbol <- rownames(gene_table)
    fgsea_results$Enriched_GENE <- .annotate_genes(
        fgsea_results, gene_table, ChemicalName_GeneEntrezIds
    )
    fgsea_results
}

#' Annotate GSEA results with enriched gene symbols
#' @return A character vector of comma-separated gene symbols per chemical.
#' @keywords internal
.annotate_genes <- function(results, gene_table, chemical_sets) {
    vapply(
        results$ChemicalID,
        function(x) {
            matched <- gene_table$EntrezID %in% chemical_sets[[x]]
            paste(gene_table[matched, "Symbol"], collapse = ", ")
        },
        character(1)
    )
}
