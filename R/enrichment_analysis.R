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
#'     Input: data frame with column \code{EntrezID} (character or
#'     numeric Entrez gene IDs) and an optional numeric value column.}
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
#'       columns, \code{EntrezID} (character or numeric Entrez gene IDs)
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
#'     of enrichment results sorted by \code{PValueAdjusted} ascending. All
#'     three methods share the leading columns \code{ChemicalID},
#'     \code{ChemicalName}, \code{Method}, \code{PValue},
#'     \code{PValueAdjusted}; method-specific extras follow (see the
#'     package vignette for the full per-method schema).
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
#'     EntrezID = c("7124", "3569", "7157", "672", "1956"),
#'     pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
#' )
#' ora_results <- enrichment_CTD(genes, method = "ORA")
#' gsea_results <- enrichment_CTD(genes, method = "GSEA")
#'
#' # CAMERA / GSVA: expression matrix + design + contrast.
#' # Uses the bundled GSE311566 subset
#' # (Dex vs DMSO, female PBMCs; see inst/extdata/README.md).
#' gse <- readRDS(system.file(
#'     "extdata", "GSE311566_subset.rds", package = "ctdR"
#' ))
#' expr <- gse$expr
#' grp  <- gse$coldata$group
#' d    <- model.matrix(~ grp)
#' camera_results <- enrichment_CTD(expr, method = "CAMERA",
#'     design = d, contrast = 2)
#'
#' # GSVA: per-sample enrichment scores
#' gsva_scores <- enrichment_CTD(expr, method = "GSVA")
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
    cache_dir <- rappdirs::user_cache_dir("ctdR")

    .validate_enrichment_args(x, method, design, contrast,
        pAdjustMethod, cache_dir)

    e <- new.env(parent = emptyenv())
    load(file = file.path(cache_dir, "chemicals.rda"), envir = e)
    chemicals <- e$chemicals

    switch(method,
        ORA = .run_ora(x, chemicals, cache_dir, pAdjustMethod),
        GSEA = .run_gsea(x, chemicals, cache_dir, pAdjustMethod),
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

#' Validate user inputs for enrichment_CTD()
#'
#' Checks pAdjustMethod, CTD cache presence, x shape per method, and
#' CAMERA-specific design/contrast requirements. Centralizing here keeps
#' enrichment_CTD() below the BiocCheck 50-line recommendation.
#'
#' @param x The user input (data.frame or matrix).
#' @param method Already normalized via match.arg().
#' @param design Design matrix for CAMERA (or NULL).
#' @param contrast Contrast spec for CAMERA (or NULL).
#' @param pAdjustMethod Multiple-testing correction name.
#' @param cache_dir Directory expected to hold CTD \code{.rda} files.
#'
#' @return Invisibly \code{TRUE} on success; otherwise \code{stop()}s
#'   with an informative message.
#' @keywords internal
.validate_enrichment_args <- function(x, method, design, contrast,
    pAdjustMethod, cache_dir) {
    valid_methods <- c("BH", "bonferroni", "fdr", "none")
    if (!pAdjustMethod %in% valid_methods) {
        stop("'pAdjustMethod' must be one of: ",
            paste(valid_methods, collapse = ", "),
            ". Got '", pAdjustMethod, "'",
            call. = FALSE
        )
    }

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
    invisible(TRUE)
}

#' Run ORA branch of enrichment analysis
#'
#' Loads the cached chemical->gene-symbol mapping, maps Entrez IDs to
#' symbols, runs ORA, and merges chemical names.
#'
#' @param x Data frame with column \code{EntrezID}.
#' @param chemicals_meta Data frame with \code{ChemicalID} and
#'   \code{ChemicalName} columns.
#' @param cache_dir Directory holding cached CTD \code{.rda} files.
#' @param pAdjustMethod Multiple-testing correction name.
#'
#' @return A data frame of ORA enrichment results.
#' @keywords internal
.run_ora <- function(x, chemicals_meta, cache_dir, pAdjustMethod) {
    e <- new.env(parent = emptyenv())
    load(
        file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"),
        envir = e
    )

    gene_symbols <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = x$EntrezID,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    )[, "SYMBOL"]

    res <- ora(
        e$ChemicalName_GeneSymbols, gene_symbols,
        pAdjustMethod = pAdjustMethod
    )

    .format_enrichment_result(res, chemicals_meta, pAdjustMethod,
        method = "ORA",
        rename = c(
            pvalue         = "PValue",
            qvalue         = "QValue",
            BgRatio        = "BackgroundRatio",
            geneID         = "EnrichedGenes",
            foldEnrichment = "FoldEnrichment"
        ),
        drop = c("p.adjust", "Description")
    )
}

#' Run GSEA branch of enrichment analysis
#'
#' Loads the cached chemical->Entrez-ID mapping, prepares the ranked
#' input, calls the \code{\link{gsea}} engine, and feeds the result
#' through the shared \code{\link{.format_enrichment_result}}.
#'
#' @param x Data frame with column \code{EntrezID} and a numeric
#'   value column used for ranking.
#' @param chemicals_meta Data frame with \code{ChemicalID} and
#'   \code{ChemicalName} columns.
#' @param cache_dir Directory holding cached CTD \code{.rda} files.
#' @param pAdjustMethod Multiple-testing correction name.
#'
#' @return A data frame of GSEA enrichment results, formatted by
#'   \code{\link{.format_enrichment_result}}.
#' @keywords internal
.run_gsea <- function(x, chemicals_meta, cache_dir, pAdjustMethod) {
    e <- new.env(parent = emptyenv())
    load(
        file = file.path(cache_dir,
            "ChemicalName_GeneEntrezIds.rda"),
        envir = e
    )
    gene_table <- as.data.frame(x)
    gene_table <- gene_table[!is.na(gene_table$EntrezID), ]
    res <- gsea(e$ChemicalName_GeneEntrezIds, gene_table)

    .format_enrichment_result(res, chemicals_meta, pAdjustMethod,
        method = "GSEA",
        rename = c(
            pval           = "PValue",
            ES             = "EnrichmentScore",
            NES            = "NormalizedEnrichmentScore",
            size           = "GeneSetSize",
            leadingEdge    = "LeadingEdge",
            foldEnrichment = "FoldEnrichment",
            Enriched_GENE  = "EnrichedGenes"
        ),
        drop = c("padj")
    )
}

#' Canonicalize and order an enrichment result data frame
#'
#' Shared post-processing for the ORA / GSEA / CAMERA runners. Each
#' runner hands in a data frame that already has the harmonized columns
#' \code{ChemicalID} and \code{pvalue}; this helper computes
#' \code{padj}, joins chemical metadata, applies a canonical
#' leading-column order, sorts by \code{padj} ascending, and drops row
#' names. Centralizing this step keeps the runners small and the public
#' output schema consistent across methods.
#'
#' @param res Data frame as returned by an enrichment engine
#'   (\code{ora()}, \code{gsea()}, \code{limma::camera()}), still
#'   carrying the engine's native column names.
#' @param chemicals_meta Data frame with \code{ChemicalID} and
#'   \code{ChemicalName} columns.
#' @param pAdjustMethod Method passed to \code{\link[stats]{p.adjust}}.
#' @param method Method label (one of \code{"ORA"}, \code{"GSEA"},
#'   \code{"CAMERA"}) stamped into the new \code{Method} column so
#'   results can be \code{rbind}'d across methods without losing
#'   provenance.
#' @param rename Named character vector mapping
#'   \code{old_engine_colname = "NewCanonicalName"}. Applied before the
#'   metadata merge. Centralizes the engine->canonical-schema mapping
#'   so each engine can keep its native column names.
#' @param drop Character vector of engine column names to remove
#'   before the merge / sort.
#'
#' @return A data frame whose leading columns are
#'   \code{ChemicalID, ChemicalName, Method, PValue, PValueAdjusted}
#'   (whichever of those are present after the rename), followed by
#'   the remaining columns in their original order, sorted by
#'   \code{PValueAdjusted} ascending.
#' @keywords internal
.format_enrichment_result <- function(res, chemicals_meta, pAdjustMethod,
    method, rename = NULL, drop = NULL) {
    if (length(drop)) {
        res <- res[, !colnames(res) %in% drop, drop = FALSE]
    }
    for (old in names(rename)) {
        colnames(res)[colnames(res) == old] <- rename[[old]]
    }

    # rep() instead of scalar assignment so the empty-result case
    # (engine returned 0 rows) doesn't trigger
    # "replacement has 1 row, data has 0".
    res$Method <- rep(method, nrow(res))
    res$PValueAdjusted <- stats::p.adjust(res$PValue, method = pAdjustMethod)
    res <- merge(res, chemicals_meta, by = "ChemicalID", all.x = TRUE)

    front <- c(
        "ChemicalID", "ChemicalName", "Method",
        "PValue", "PValueAdjusted"
    )
    front <- intersect(front, colnames(res))
    other <- setdiff(colnames(res), front)
    res <- res[, c(front, other), drop = FALSE]
    res <- res[order(res$PValueAdjusted), , drop = FALSE]
    rownames(res) <- NULL
    res
}
