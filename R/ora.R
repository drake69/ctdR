#' @title Over-Representation Analysis (ORA)
#'
#' @description
#' Internal function that performs Over-Representation Analysis using
#' \code{\link[clusterProfiler]{enricher}}. Tests whether the input gene list
#' is significantly enriched for genes associated with each CTD chemical.
#'
#' @param ChemicalName_GeneSymbols A data frame with two columns
#'   (\code{term}, \code{gene}) mapping CTD chemical IDs to HGNC gene symbols.
#'   This serves as the TERM2GENE input for \code{clusterProfiler::enricher}.
#' @param entrez_ids Character vector of HGNC gene symbols to test for
#'   enrichment.
#'
#' @return A data frame with columns: \code{ID}, \code{GeneRatio},
#'   \code{BgRatio}, \code{pvalue}, \code{padj}, \code{qvalue}, \code{geneID},
#'   \code{Count}, and \code{foldEnrichment}. Returns an empty data frame with
#'   the same structure if no enrichment is found.
#'
#' @keywords internal
ora <- function(ChemicalName_GeneSymbols, entrez_ids) {

  gene_list <- entrez_ids
  gene2pathway <- ChemicalName_GeneSymbols

  # Utilizza l'enricher per l'analisi ORA
  ora_results <- clusterProfiler::enricher(gene = gene_list, TERM2GENE = gene2pathway)

  if (is.null(ora_results)) {
    return(data.frame(ID = character(), GeneRatio = character(), BgRatio = character(),
                      pvalue = numeric(), padj = numeric(), qvalue = numeric(),
                      geneID = character(), Count = integer(), foldEnrichment = numeric()))
  }

  ora_results <- ora_results@result

  # rename p.adjust as padj
  colnames(ora_results)[colnames(ora_results)=="p.adjust"] <- "padj"

  # calculate foldEnrichment using BgRatio and GeneRatio
  # dplyr::mutate(ora_results, foldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
  ora_results$foldEnrichment <- DOSE::parse_ratio(ora_results$GeneRatio) / DOSE::parse_ratio(ora_results$BgRatio)

  return(ora_results)
}
