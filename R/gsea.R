#' @title Gene Set Enrichment Analysis (GSEA)
#'
#' @description
#' Internal function that performs Gene Set Enrichment Analysis using
#' \code{\link[fgsea]{fgsea}}. Detects chemicals whose known target genes
#' cluster toward the extremes of a ranked gene list.
#'
#' @param ChemicalName_GeneSymbols A named list where each element is a
#'   character vector of Entrez gene IDs associated with a CTD chemical. Names
#'   are CTD chemical IDs.
#' @param entrez_ids A data frame with at least two columns:
#'   \describe{
#'     \item{\code{entrez_ids}}{Entrez gene IDs (character).}
#'     \item{(second column)}{Numeric values used for ranking (e.g. p-values).
#'       Converted internally to \code{-log10(value)}
#'       for the ranking statistic.}
#'   }
#' @param chemicals A data frame with columns \code{ChemicalID} and
#'   \code{ChemicalName}, used to annotate results with human-readable names.
#' @param pAdjustMethod Character. Method for multiple testing correction
#'   (default \code{"BH"}). Passed to \code{\link[stats]{p.adjust}} to
#'   recalculate the \code{padj} column returned by \code{fgsea}.
#'
#' @return A data frame with columns: \code{ChemicalID}, \code{pval},
#'   \code{padj}, \code{ES}, \code{NES}, \code{size}, \code{leadingEdge},
#'   \code{ChemicalName}, \code{foldEnrichment}, and \code{Enriched_GENE}.
#'   Results are sorted by \code{padj} in ascending order.
#'
#' @keywords internal
gsea <- function(ChemicalName_GeneSymbols, entrez_ids,
    chemicals, pAdjustMethod = "BH") {
    stats <- -log10(entrez_ids[, 2] + 1e-10)
    names(stats) <- entrez_ids$entrez_ids
    stats <- sort(stats, decreasing = TRUE)

    fgsea_results <- fgsea::fgsea(
        pathways = ChemicalName_GeneSymbols, stats = stats
    )
    fgsea_results <- as.data.frame(fgsea_results)
    fgsea_results$padj <- stats::p.adjust(
        fgsea_results$pval, method = pAdjustMethod
    )
    fgsea_results <- fgsea_results[order(fgsea_results$padj), ]
    fgsea_results$foldEnrichment <- (
        abs(fgsea_results$ES) / mean(fgsea_results$ES)
    )
    names(fgsea_results)[1] <- "ChemicalID"
    fgsea_results <- merge(
        fgsea_results, chemicals, by = "ChemicalID"
    )
    entrez_ids$Symbol <- rownames(entrez_ids)
    fgsea_results$Enriched_GENE <- .annotate_genes(
        fgsea_results, entrez_ids, ChemicalName_GeneSymbols
    )
    fgsea_results
}

#' Annotate GSEA results with enriched gene symbols
#' @return A character vector of comma-separated gene symbols per chemical.
#' @keywords internal
.annotate_genes <- function(results, entrez_ids, pathways) {
    vapply(
        results$ChemicalID,
        function(x) {
            matched <- entrez_ids$entrez_ids %in% pathways[[x]]
            paste(entrez_ids[matched, "Symbol"], collapse = ", ")
        },
        character(1)
    )
}
