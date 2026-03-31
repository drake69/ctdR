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
#'     larger than expected by chance.
#'     Uses \code{\link[clusterProfiler]{enricher}}.
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
#'   \code{ChemicalID} \tab CTD chemical identifier \cr
#'   \code{ChemicalName}   \tab Human-readable chemical name \cr
#'   \code{GeneRatio} \tab Proportion of input genes in set \cr
#'   \code{BgRatio}        \tab Background ratio \cr
#'   \code{pvalue}         \tab Raw p-value \cr
#'   \code{padj} \tab Adjusted p-value \cr
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
#' @param pAdjustMethod Character. Method for multiple testing correction.
#'   One of \code{"BH"} (Benjamini-Hochberg, default), \code{"bonferroni"},
#'   \code{"fdr"} (alias for BH), or \code{"none"}.
#'   Passed to \code{\link[clusterProfiler]{enricher}} for ORA and to
#'   \code{\link[stats]{p.adjust}} for GSEA.
#'
#' @return A data frame of enrichment results (see Details for columns).
#'
#' @seealso \code{\link{import_CTD}} to import and cache the CTD data before
#'   first use.
#'
#' @examples
#' # Import the bundled sample data first:
#' sample_file <- system.file(
#'     "extdata", "CTD_chem_gene_ixns_sample.csv",
#'     package = "ctdR"
#' )
#' import_CTD(sample_file)
#'
#' # Prepare a gene list (Entrez IDs with p-values):
#' genes <- data.frame(
#'     entrez_ids = c("7124", "3569", "7157", "672", "1956"),
#'     pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
#' )
#'
#' # Over-Representation Analysis
#' ora_results <- enrichment_CTD(genes, method = "ORA")
#' head(ora_results)
#'
#' # Gene Set Enrichment Analysis
#' gsea_results <- enrichment_CTD(genes, method = "GSEA")
#' head(gsea_results)
#'
#' @importFrom rappdirs user_cache_dir
#' @export
enrichment_CTD <- function(entrez_ids, method = "ORA", pAdjustMethod = "BH") {
    valid_methods <- c("BH", "bonferroni", "fdr", "none")
    if (!pAdjustMethod %in% valid_methods) {
        stop("'pAdjustMethod' must be one of: ",
            paste(valid_methods, collapse = ", "),
            ". Got '", pAdjustMethod, "'",
            call. = FALSE
        )
    }

    cache_dir <- rappdirs::user_cache_dir("ctdR")

    rda_path <- file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda")
    if (!file.exists(rda_path)) {
        stop(
            "CTD data not found. Please:\n",
            "  1. Download CTD_chem_gene_ixns.csv.gz from\n",
            "     https://ctdbase.org/reports/",
            "CTD_chem_gene_ixns.csv.gz\n",
            "  2. Decompress: gunzip CTD_chem_gene_ixns.csv.gz\n",
            "  3. Run import_CTD(\"path/to/file.csv\")",
            call. = FALSE
        )
    }

    load(file = file.path(cache_dir, "chemicals.rda"))

    if (method == "ORA") {
        load(file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"))
        fgsea_results <- .run_ora(
            entrez_ids, ChemicalName_GeneSymbols,
            chemicals, pAdjustMethod
        )
    } else {
        load(file = file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"))
        entrez_ids <- as.data.frame(entrez_ids)
        entrez_ids <- entrez_ids[!is.na(entrez_ids$entrez_ids), ]
        fgsea_results <- gsea(
            ChemicalName_GeneEntrezIds, entrez_ids,
            chemicals,
            pAdjustMethod = pAdjustMethod
        )
    }

    return(fgsea_results)
}

#' Run ORA branch of enrichment analysis
#' @return A data frame of ORA enrichment results.
#' @keywords internal
.run_ora <- function(entrez_ids, ChemicalName_GeneSymbols,
    chemicals, pAdjustMethod) {
    gene_symbols <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = entrez_ids$entrez_ids,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    )
    gene_symbols <- gene_symbols[, "SYMBOL"]
    results <- ora(
        ChemicalName_GeneSymbols, gene_symbols,
        pAdjustMethod = pAdjustMethod
    )
    results <- merge(
        results, chemicals,
        by.y = "ChemicalID", by.x = "ID"
    )
    id_col <- colnames(results) == "ID"
    colnames(results)[id_col] <- "ChemicalID"
    keep <- !colnames(results) %in% "Description"
    results[, keep]
}
