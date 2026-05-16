#' @title Chemical Enrichment Analysis Using CTD
#'
#' @description
#' Identifies chemicals whose known gene targets are significantly enriched in
#' a user-supplied gene list or expression matrix, using data from the
#' Comparative Toxicogenomics Database (CTD).
#'
#' Four methods are available:
#' \describe{
#'   \item{\strong{ORA}}{Over-Representation Analysis (default). Tests whether
#'     the overlap between your gene list and each chemical's target genes is
#'     larger than expected by chance.
#'     Uses \code{\link[clusterProfiler]{enricher}}.
#'     Input: data frame of Entrez IDs (column \code{entrez_ids}) with an
#'     optional numeric value column.}
#'   \item{\strong{GSEA}}{Gene Set Enrichment Analysis. Uses a ranked gene list
#'     (ranked by the numeric column in the input, e.g. p-values or fold
#'     changes) to detect chemicals whose targets cluster toward the top or
#'     bottom of the ranking. Uses \code{\link[fgsea]{fgsea}}.}
#'   \item{\strong{CAMERA}}{Competitive gene-set test accounting for
#'     inter-gene correlation. Uses \code{\link[limma]{camera}}.
#'     Input: a numeric expression matrix (genes x samples) plus a design
#'     matrix and a contrast.}
#'   \item{\strong{GSVA}}{Gene Set Variation Analysis (sample-level scoring).
#'     Returns a matrix of per-sample enrichment scores for each chemical,
#'     suitable for downstream clustering or association testing.
#'     Uses \code{\link[GSVA]{gsva}}.
#'     Input: a numeric expression matrix (genes x samples).}
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
#' @param x The input. Its expected type depends on \code{method}:
#'   \itemize{
#'     \item For \code{"ORA"} and \code{"GSEA"}: a data frame with at least two
#'       columns, \code{entrez_ids} (character or numeric Entrez gene IDs)
#'       and a numeric value column (e.g. p-value, log fold-change). The
#'       numeric column is used for ranking in GSEA and ignored in ORA.
#'     \item For \code{"CAMERA"} and \code{"GSVA"}: a numeric expression
#'       matrix with genes in rows and samples in columns. \code{rownames(x)}
#'       must be either Entrez IDs or HGNC SYMBOLs.
#'   }
#' @param method Character. Enrichment method: \code{"ORA"} (default),
#'   \code{"GSEA"}, \code{"CAMERA"}, or \code{"GSVA"}.
#' @param design Design matrix (required when \code{method = "CAMERA"}).
#' @param contrast Contrast specification for \code{\link[limma]{camera}}
#'   (column number, column name, or numeric vector). Required when
#'   \code{method = "CAMERA"}.
#' @param id_type Either \code{"entrez"}, \code{"symbol"}, or \code{NULL}
#'   (default) for auto-detection from \code{rownames(x)}. Only used when
#'   \code{method} is \code{"CAMERA"} or \code{"GSVA"}.
#' @param pAdjustMethod Character. Method for multiple testing correction:
#'   one of \code{"BH"} (Benjamini-Hochberg, default), \code{"bonferroni"},
#'   \code{"fdr"} (alias for BH), or \code{"none"}. Not used for
#'   \code{method = "GSVA"} (which returns scores rather than p-values).
#' @param ... Additional arguments forwarded to the underlying engine:
#'   \code{\link[clusterProfiler]{enricher}} for ORA,
#'   \code{\link[fgsea]{fgsea}} for GSEA,
#'   \code{\link[limma]{camera}} for CAMERA,
#'   \code{\link[GSVA]{gsva}} for GSVA.
#'
#' @return
#' \itemize{
#'   \item For \code{"ORA"}, \code{"GSEA"}, and \code{"CAMERA"}: a data frame
#'     of enrichment results sorted by adjusted p-value (\code{padj}).
#'   \item For \code{"GSVA"}: a numeric matrix of enrichment scores with
#'     chemicals (CTD chemical IDs) in rows and samples in columns.
#' }
#'
#' @seealso \code{\link{import_CTD}} to import and cache the CTD data;
#'   \code{\link{plot_CTD}} to visualize results.
#'
#' @examples
#' # Import the bundled sample data first:
#' sample_file <- system.file(
#'     "extdata", "CTD_chem_gene_ixns_sample.csv",
#'     package = "ctdR"
#' )
#' import_CTD(sample_file)
#'
#' # ORA / GSEA: prepare a gene list with Entrez IDs and a numeric value
#' genes <- data.frame(
#'     entrez_ids = c("7124", "3569", "7157", "672", "1956"),
#'     pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
#' )
#' ora_results <- enrichment_CTD(genes, method = "ORA")
#' gsea_results <- enrichment_CTD(genes, method = "GSEA")
#'
#' # CAMERA: expression matrix + design + contrast
#' set.seed(1)
#' expr <- matrix(rnorm(20 * 6), nrow = 20,
#'     dimnames = list(c("7124","3569","7157","672","1956",
#'         "1017","1019","207","208","595",
#'         "894","983","1029","1869","5599",
#'         "5290","4609","6347","3458","2099"),
#'         paste0("S", 1:6)))
#' grp <- factor(rep(c("ctrl","treat"), each = 3))
#' d <- model.matrix(~ grp)
#' \donttest{
#' camera_results <- enrichment_CTD(expr, method = "CAMERA",
#'     design = d, contrast = 2)
#'
#' # GSVA: per-sample enrichment scores
#' gsva_scores <- enrichment_CTD(expr, method = "GSVA")
#' }
#'
#' @importFrom rappdirs user_cache_dir
#' @export
enrichment_CTD <- function(x,
    method = c("ORA", "GSEA", "CAMERA", "GSVA"),
    design = NULL,
    contrast = NULL,
    id_type = NULL,
    pAdjustMethod = "BH",
    ...) {
    if (missing(x)) {
        stop("Argument 'x' is required.", call. = FALSE)
    }
    method <- match.arg(method)

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

    if (method %in% c("ORA", "GSEA")) {
        if (!is.data.frame(x)) {
            stop("For method = '", method,
                "', 'x' must be a data.frame with Entrez IDs ",
                "and a numeric value column.",
                call. = FALSE
            )
        }
    } else {
        if (!is.matrix(x) || !is.numeric(x)) {
            stop("For method = '", method,
                "', 'x' must be a numeric matrix (genes x samples).",
                call. = FALSE
            )
        }
        if (method == "CAMERA" &&
                (is.null(design) || is.null(contrast))) {
            stop(
                "For method = 'CAMERA', 'design' and 'contrast' are required.",
                call. = FALSE
            )
        }
    }

    e <- new.env(parent = emptyenv())
    load(file = file.path(cache_dir, "chemicals.rda"), envir = e)
    chemicals <- e$chemicals

    switch(method,
        ORA = {
            e2 <- new.env(parent = emptyenv())
            load(
                file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"),
                envir = e2
            )
            .run_ora(x, e2$ChemicalName_GeneSymbols,
                chemicals, pAdjustMethod)
        },
        GSEA = {
            e2 <- new.env(parent = emptyenv())
            load(
                file = file.path(cache_dir,
                    "ChemicalName_GeneEntrezIds.rda"),
                envir = e2
            )
            entrez <- as.data.frame(x)
            entrez <- entrez[!is.na(entrez$entrez_ids), ]
            gsea(e2$ChemicalName_GeneEntrezIds, entrez,
                chemicals, pAdjustMethod = pAdjustMethod)
        },
        CAMERA = .run_camera(
            expr = x, design = design, contrast = contrast,
            id_type = id_type, pAdjustMethod = pAdjustMethod,
            chemicals_meta = chemicals,
            cache_dir = cache_dir,
            ...
        ),
        GSVA = .run_gsva(
            expr = x,
            id_type = id_type,
            cache_dir = cache_dir,
            ...
        )
    )
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
