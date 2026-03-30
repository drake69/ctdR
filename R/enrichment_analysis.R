#' @title Chemical Enrichment Analysis Using CTD
#'
#' @description
#' Identifies chemicals whose known gene targets are significantly enriched in a
#' user-supplied gene list, using data from the Comparative Toxicogenomics
#' Database (CTD).
#'
#' Two methods are available:
#' \describe{
#'   \item{\strong{ORA}}{Over-Representation Analysis (default). Tests whether
#'     the overlap between your gene list and each chemical's target genes is
#'     larger than expected by chance. Uses \code{\link[clusterProfiler]{enricher}}.
#'     Input genes are mapped from Entrez IDs to HGNC symbols internally.}
#'   \item{\strong{GSEA}}{Gene Set Enrichment Analysis. Uses a ranked gene list
#'     (ranked by the numeric column in \code{entrez_ids}, e.g. p-values or fold
#'     changes) to detect chemicals whose targets cluster toward the top or
#'     bottom of the ranking. Uses \code{\link[fgsea]{fgsea}}.}
#' }
#'
#' @details
#' Before calling this function you must import the CTD data once with
#' \code{\link{import_CTD}}. If the cached data is not found, the function
#' stops with an informative error message.
#'
#' @section Data Licensing Disclaimer:
#' This package does \strong{not} bundle or redistribute any CTD data. The
#' Comparative Toxicogenomics Database is maintained by NC State University and
#' its data are subject to specific licensing terms. Users must download the
#' data directly from \url{https://ctdbase.org} and comply with the CTD Terms
#' of Service (\url{https://ctdbase.org/about/legal.jsp}).
#'
#' \strong{ORA output columns:}
#' \tabular{ll}{
#'   \code{ChemicalID}     \tab CTD chemical identifier (e.g. \code{"D000082"}) \cr
#'   \code{ChemicalName}   \tab Human-readable chemical name \cr
#'   \code{GeneRatio}      \tab Proportion of input genes in the chemical's set \cr
#'   \code{BgRatio}        \tab Background ratio \cr
#'   \code{pvalue}         \tab Raw p-value \cr
#'   \code{padj}           \tab Adjusted p-value (BH method) \cr
#'   \code{foldEnrichment} \tab GeneRatio / BgRatio \cr
#'   \code{geneID}         \tab Enriched gene symbols \cr
#'   \code{Count}          \tab Number of overlapping genes \cr
#' }
#'
#' \strong{GSEA output columns:}
#' \tabular{ll}{
#'   \code{ChemicalID}     \tab CTD chemical identifier \cr
#'   \code{ChemicalName}   \tab Human-readable chemical name \cr
#'   \code{pval}           \tab Raw p-value \cr
#'   \code{padj}           \tab Adjusted p-value \cr
#'   \code{ES}             \tab Enrichment score \cr
#'   \code{NES}            \tab Normalized enrichment score \cr
#'   \code{size}           \tab Size of the gene set \cr
#'   \code{leadingEdge}    \tab Leading-edge gene subset \cr
#'   \code{foldEnrichment} \tab |ES| / mean(ES) \cr
#'   \code{Enriched_GENE}  \tab Comma-separated enriched gene symbols \cr
#' }
#'
#' @param entrez_ids A data frame with at least two columns:
#'   \describe{
#'     \item{\code{entrez_ids}}{Character or numeric Entrez gene IDs.}
#'     \item{(second column)}{A numeric value associated with each gene (e.g.
#'       p-value, log fold-change). Used for ranking in GSEA; ignored in ORA.}
#'   }
#' @param method Character. Enrichment method to use: \code{"ORA"} (default) for
#'   Over-Representation Analysis or \code{"GSEA"} for Gene Set Enrichment
#'   Analysis.
#'
#' @return A data frame of enrichment results (see Details for columns).
#'
#' @seealso \code{\link{import_CTD}} to import and cache the CTD data before
#'   first use.
#'
#' @examples
#' # enrichment_CTD checks for cached CTD data before running:
#' genes <- data.frame(
#'   entrez_ids = c("7124", "3569", "7157", "672", "1956"),
#'   pvalue     = c(0.001, 0.003, 0.01, 0.02, 0.05)
#' )
#' tryCatch(
#'   enrichment_CTD(genes, method = "ORA"),
#'   error = function(e) message(e$message)
#' )
#'
#' \donttest{
#' # Full workflow (requires imported CTD data):
#' # import_CTD("~/Downloads/CTD_chem_gene_ixns.csv")
#'
#' # Over-Representation Analysis
#' ora_results <- enrichment_CTD(genes, method = "ORA")
#' head(ora_results)
#'
#' # Gene Set Enrichment Analysis
#' gsea_results <- enrichment_CTD(genes, method = "GSEA")
#' head(gsea_results)
#' }
#'
#' @importFrom rappdirs user_cache_dir
#' @export
enrichment_CTD <- function(entrez_ids, method="ORA")
{
  cache_dir <- rappdirs::user_cache_dir("ctdR")

  if (!file.exists(file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"))) {
    stop("CTD data not found. Please:\n",
         "  1. Download CTD_chem_gene_ixns.csv.gz from https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz\n",
         "  2. Decompress it: gunzip CTD_chem_gene_ixns.csv.gz\n",
         "  3. Run import_CTD(\"path/to/CTD_chem_gene_ixns.csv\") once to import and cache the data",
         call. = FALSE)
  }

  load(file = file.path(cache_dir, "chemicals.rda"))

  if(method=="ORA")
  {
    load(file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"))
    gene_symbols <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = entrez_ids$entrez_ids,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first")
    gene_symbols <- gene_symbols[,"SYMBOL"]
    fgsea_results <- ora(ChemicalName_GeneSymbols,gene_symbols)
    # Add ChemicalName
    fgsea_results <- merge(fgsea_results, chemicals, by.y="ChemicalID", by.x="ID")
    # rename ID as ChemicalID
    colnames(fgsea_results)[colnames(fgsea_results)=="ID"] <- "ChemicalID"
    # remove Description column
    fgsea_results <- fgsea_results[,!colnames(fgsea_results) %in% c("Description")]
  }
  else
  {
    load(file = file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"))
    entrez_ids <- as.data.frame(entrez_ids)
    # Remove NA values
    entrez_ids <- entrez_ids[!is.na(entrez_ids$entrez_ids),]
    fgsea_results <- gsea(ChemicalName_GeneEntrezIds, entrez_ids, chemicals)
  }

  return(fgsea_results)
}
