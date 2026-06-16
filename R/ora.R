#' @title Over-Representation Analysis (ORA)
#'
#' @description
#' Internal engine that performs Over-Representation Analysis via
#' \code{\link[clusterProfiler]{enricher}}. Returns the raw enricher
#' result table with one ctdR-added column (\code{foldEnrichment}).
#'
#' Column renaming, multiple-testing correction, chemical-name join,
#' canonical ordering and sort are applied downstream by
#' \code{\link{.format_enrichment_result}} so the engine stays close
#' to its underlying tool's native vocabulary.
#'
#' @param ChemicalName_GeneSymbols A data frame with two columns
#'   (\code{term}, \code{gene}) mapping CTD chemical IDs to HGNC gene
#'   symbols. Serves as the TERM2GENE input for
#'   \code{clusterProfiler::enricher}.
#' @param gene_symbols Character vector of HGNC gene symbols to test
#'   for enrichment.
#' @param pAdjustMethod Character. Method for multiple testing
#'   correction (default \code{"BH"}). Passed to
#'   \code{\link[clusterProfiler]{enricher}}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[clusterProfiler]{enricher}}, e.g.:
#'   \describe{
#'     \item{\code{universe}}{Character vector of background gene symbols
#'       (default: all genes in the TERM2GENE table). Set to
#'       \code{rownames(expr)} or the full tested gene list to restrict
#'       the background to measured genes only.}
#'     \item{\code{minGSSize}}{Minimum gene set size after intersection
#'       with the universe (default 1).}
#'     \item{\code{maxGSSize}}{Maximum gene set size (default 500).}
#'   }
#'
#' @return A data frame with \code{clusterProfiler::enricher}'s native
#'   columns (\code{ID}, \code{GeneRatio}, \code{BgRatio},
#'   \code{pvalue}, \code{p.adjust}, \code{qvalue}, \code{geneID},
#'   \code{Count}, \code{Description}) plus a ctdR-added
#'   \code{foldEnrichment} column. Returns an empty data frame with
#'   the same structure if no enrichment is found.
#'
#' @keywords internal
ora <- function(ChemicalName_GeneSymbols, gene_symbols,
    pAdjustMethod = "BH", ...) {
    empty <- data.frame(
        ChemicalID = character(), Description = character(),
        GeneRatio = character(), BgRatio = character(),
        pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
        geneID = character(), Count = integer(),
        foldEnrichment = numeric()
    )

    ora_results <- clusterProfiler::enricher(
        gene = gene_symbols,
        TERM2GENE = ChemicalName_GeneSymbols,
        pAdjustMethod = pAdjustMethod,
        ...
    )

    if (is.null(ora_results)) return(empty)

    res <- ora_results@result
    # clusterProfiler labels its primary key "ID", but in ctdR those
    # values are CTD chemical IDs. Lift to the semantically correct
    # name here so callers never see a generic "ID".
    colnames(res)[colnames(res) == "ID"] <- "ChemicalID"
    res$foldEnrichment <-
        .parse_ratio(res$GeneRatio) / .parse_ratio(res$BgRatio)
    res
}

#' Parse "n/d" ratio strings into numeric values
#' @param x Character vector of strings of the form "n/d".
#' @return Numeric vector of n/d ratios.
#' @keywords internal
.parse_ratio <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]),
        numeric(1))
}
