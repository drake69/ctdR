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
#'       Converted internally to \code{-log10(value)} for the ranking statistic.}
#'   }
#' @param chemicals A data frame with columns \code{ChemicalID} and
#'   \code{ChemicalName}, used to annotate results with human-readable names.
#'
#' @return A data frame with columns: \code{ChemicalID}, \code{pval},
#'   \code{padj}, \code{ES}, \code{NES}, \code{size}, \code{leadingEdge},
#'   \code{ChemicalName}, \code{foldEnrichment}, and \code{Enriched_GENE}.
#'   Results are sorted by \code{padj} in ascending order.
#'
#' @keywords internal
gsea <- function(ChemicalName_GeneSymbols, entrez_ids, chemicals)
{

  # Create a named numeric vector for fgsea
  stats <- -log10(entrez_ids[,2] + 1e-10)
  names(stats) <- entrez_ids$entrez_ids
  stats <- sort(stats, decreasing = TRUE)

  # Run fgsea
  fgsea_results <- fgsea::fgsea(pathways = ChemicalName_GeneSymbols, stats = stats)

  fgsea_results <- as.data.frame(fgsea_results)

  # sort by padj asc
  fgsea_results <- fgsea_results[order(fgsea_results$padj),]

  fgsea_results$foldEnrichment <- abs(fgsea_results$ES) / mean(fgsea_results$ES)

  # change pathway column name to substance
  names(fgsea_results)[1] <- "ChemicalID"

  # Add ChemicalName
  fgsea_results <- merge(fgsea_results, chemicals, by="ChemicalID")

  entrez_ids$Symbol <- rownames(entrez_ids)
  # Add Enriched GENE
  fgsea_results$Enriched_GENE <- vapply(
    fgsea_results$ChemicalID,
    function(x) {
      paste(
        entrez_ids[
          entrez_ids$entrez_ids %in% ChemicalName_GeneSymbols[[x]],
          "Symbol"
        ],
        collapse = ", "
      )
    },
    character(1)
  )

  return(fgsea_results)

}
