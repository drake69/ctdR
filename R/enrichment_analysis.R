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
#'     larger than expected by chance. Uses a hypergeometric test
#'     computed directly on \code{\link[stats]{phyper}}.
#'     Input: data frame with column \code{EntrezID} (character or
#'     numeric Entrez gene IDs) and an optional numeric value column.}
#'   \item{\strong{GSEA}}{Gene Set Enrichment Analysis. Uses a ranked gene list
#'     (ranked by the numeric column in the input, e.g. p-values or fold
#'     changes) to detect chemicals whose targets cluster toward the top or
#'     bottom of the ranking. Uses \code{\link[fgsea]{fgsea}}.}
#'   \item{\strong{CAMERA}}{Competitive gene-set test accounting for
#'     inter-gene correlation. Uses \code{\link[limma]{camera}}.
#'     Input: a numeric expression matrix (genes x samples) or a
#'     \code{\link[SummarizedExperiment]{SummarizedExperiment}}, plus a
#'     design matrix and a contrast.}
#'   \item{\strong{GSVA}}{Gene Set Variation Analysis (sample-level scoring).
#'     Returns per-sample enrichment scores for each chemical, suitable for
#'     downstream clustering or association testing.
#'     Uses \code{\link[GSVA]{gsva}}.
#'     Input: a numeric expression matrix (genes x samples) or a
#'     \code{\link[SummarizedExperiment]{SummarizedExperiment}}, which is
#'     also the container the scores are returned in.}
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
#'       and a numeric value column (e.g. p-value). For GSEA, an optional
#'       column named \code{stat} can be added with a signed ranking statistic
#'       (e.g. the moderated t-statistic from \code{limma::eBayes()});
#'       when present it is used directly for ranking, preserving directionality
#'       and avoiding ties. When absent, the second column is transformed via
#'       \code{-log10()} with a warning. The second column is ignored by ORA.
#'     \item For \code{"CAMERA"} and \code{"GSVA"}: either a numeric
#'       expression matrix with genes in rows and samples in columns, or a
#'       \code{\link[SummarizedExperiment]{SummarizedExperiment}} whose
#'       assay holds that matrix (see \code{assay} to choose which one).
#'       \code{rownames(x)} must be either Entrez IDs or HGNC SYMBOLs.
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
#' @param pAdjustMethod Character. Multiple-testing correction applied to the
#'   raw p-values, passed to \code{\link[stats]{p.adjust}}. Any value in
#'   \code{stats::p.adjust.methods} is accepted:
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BH"} (Benjamini-Hochberg, the default), \code{"BY"},
#'   \code{"fdr"} (alias for \code{"BH"}), or \code{"none"}. See
#'   \code{\link[stats]{p.adjust}} for the meaning of each method. Not used for
#'   \code{method = "GSVA"} (which returns per-sample scores, not p-values).
#' @param gene_id_type Character. Identifier type used in the \code{EnrichedGenes}
#'   output column: \code{"symbol"} (default) returns HGNC gene symbols with
#'   Entrez ID as fallback for unmapped genes; \code{"entrez"} skips the
#'   symbol lookup and returns Entrez IDs directly. Only applies to \code{"ORA"}
#'   and \code{"GSEA"}; ignored by \code{"CAMERA"} and \code{"GSVA"}.
#' @param interaction_types Character vector of CTD \code{InteractionActions}
#'   values to retain when building gene sets, or \code{NULL} (default) to
#'   use all cached interactions. Values follow the \code{verb\^{}noun}
#'   convention used by CTD, e.g. \code{"increases\^{}expression"},
#'   \code{"decreases\^{}expression"}, \code{"affects\^{}binding"}.
#'   A gene is included in a chemical's set if \emph{any} of its recorded
#'   interaction actions matches one of the specified types. Requires that
#'   \code{import_CTD()} has been run (the filter is applied to the cached
#'   \code{ctd_interactions.rda} file). Restricting to expression interactions
#'   is recommended for RNA-seq analyses to improve biological specificity.
#' @param assay Which assay to use when \code{x} is a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}: an assay name,
#'   a positive index, or \code{NULL} (default) for the first assay. Ignored
#'   when \code{x} is a matrix or a data frame.
#' @param ... Additional arguments forwarded to the underlying engine:
#'   \code{\link{ora}} for ORA (\code{universe}, \code{minGSSize},
#'   \code{maxGSSize}; note that \code{minGSSize} defaults to 2, a value
#'   chosen for CTD rather than inherited from a general-purpose tool),
#'   \code{\link[fgsea]{fgseaMultilevel}} for GSEA (e.g. \code{minSize},
#'   \code{maxSize}, \code{nproc}),
#'   \code{\link[limma]{camera}} for CAMERA,
#'   \code{\link[GSVA]{gsva}} for GSVA (e.g. \code{minSize}, \code{maxSize}).
#'
#' @return
#' \itemize{
#'   \item For \code{"ORA"}, \code{"GSEA"}, and \code{"CAMERA"}: a data frame
#'     of enrichment results sorted by \code{PValueAdjusted} ascending. All
#'     three methods share the leading columns \code{ChemicalID},
#'     \code{ChemicalName}, \code{Method}, \code{PValue},
#'     \code{PValueAdjusted}; method-specific extras follow (see the
#'     package vignette for the full per-method schema).
#'   \item For \code{"GSVA"}: enrichment scores with chemicals (CTD chemical
#'     IDs) in rows and samples in columns. The container follows the input:
#'     a matrix in returns a numeric matrix, while a
#'     \code{\link[SummarizedExperiment]{SummarizedExperiment}} in returns a
#'     \code{SummarizedExperiment} whose assay holds the scores and whose
#'     \code{colData} is carried over from the input, so sample annotation
#'     stays attached to the results.
#' }
#'
#' @seealso \code{\link{import_CTD}} to import and cache the CTD data;
#'   \code{\link{plot_CTD}} to visualize results.
#'
#' @examples
#' # Import the bundled sample data first:
#' # Examples write to a temporary cache, so running them cannot
#' # disturb CTD data you have already imported. Set the same option
#' # yourself to keep an analysis isolated from your main cache.
#' options(ctdR.cache = tempfile())
#'
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
#' # CAMERA / GSVA: expression data plus, for CAMERA, design + contrast.
#' # Uses the bundled GSE311566 subset, a SummarizedExperiment
#' # (Dex vs DMSO, female PBMCs; see inst/extdata/README.md).
#' se <- readRDS(system.file(
#'     "extdata", "GSE311566_subset.rds", package = "ctdR"
#' ))
#' d <- model.matrix(~ se$group)
#' camera_results <- enrichment_CTD(se, method = "CAMERA",
#'     design = d, contrast = 2)
#'
#' # GSVA: SummarizedExperiment in, SummarizedExperiment out,
#' # so the sample annotation stays attached to the scores.
#' gsva_scores <- enrichment_CTD(se, method = "GSVA")
#' SummarizedExperiment::colData(gsva_scores)$group
#'
#' @export
enrichment_CTD <- function(x,
    method = c("ORA", "GSEA", "CAMERA", "GSVA"),
    design = NULL,
    contrast = NULL,
    id_type = NULL,
    pAdjustMethod = "BH",
    interaction_types = NULL,
    gene_id_type = c("symbol", "entrez"),
    assay = NULL,
    ...) {
    gene_id_type <- match.arg(gene_id_type)
    if (missing(x)) {
        stop("Argument 'x' is required.", call. = FALSE)
    }
    method <- match.arg(method)
    cache_dir <- .ctd_cache_dir()

    .validate_enrichment_args(x, method, design, contrast,
        pAdjustMethod, cache_dir)

    bfc <- .ctd_bfc(cache_dir)
    chemicals <- .ctd_cache_load(bfc, "chemicals")
    # Read once here rather than in each runner: every method carries the
    # same record, and the runners differ only in the container it goes on.
    provenance <- .ctd_provenance_cached(bfc)

    res <- switch(method,
        ORA = .run_ora(x, chemicals, cache_dir, pAdjustMethod,
                       interaction_types = interaction_types,
                       gene_id_type = gene_id_type, ...),
        GSEA = .run_gsea(x, chemicals, cache_dir, pAdjustMethod,
                         interaction_types = interaction_types,
                         gene_id_type = gene_id_type, ...),
        CAMERA = .run_camera(
            expr = x, design = design, contrast = contrast,
            id_type = id_type, pAdjustMethod = pAdjustMethod,
            chemicals_meta = chemicals,
            cache_dir = cache_dir,
            interaction_types = interaction_types,
            assay = assay,
            ...
        ),
        GSVA = .run_gsva(
            expr = x,
            id_type = id_type,
            cache_dir = cache_dir,
            interaction_types = interaction_types,
            assay = assay,
            ...
        )
    )
    .attach_provenance(res, provenance)
}

#' Validate user inputs for enrichment_CTD()
#'
#' Checks pAdjustMethod, CTD cache presence, x shape per method, and
#' CAMERA-specific design/contrast requirements. Centralizing here keeps
#' enrichment_CTD() below the BiocCheck 50-line recommendation.
#'
#' @param x The user input (data.frame, matrix, or SummarizedExperiment).
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
    if (!pAdjustMethod %in% stats::p.adjust.methods) {
        stop("'pAdjustMethod' must be one of: ",
            paste(stats::p.adjust.methods, collapse = ", "),
            ". Got '", pAdjustMethod, "'",
            call. = FALSE
        )
    }

    if (!.ctd_cache_has(.ctd_bfc(cache_dir), "ChemicalName_GeneEntrezIds")) {
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
        if (!.is_se(x) && !(is.matrix(x) && is.numeric(x))) {
            stop("For method = '", method,
                "', 'x' must be a numeric matrix (genes x samples) ",
                "or a SummarizedExperiment.",
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
#' @param interaction_types Character vector of CTD \code{InteractionActions}
#'   values to retain when building gene sets, or \code{NULL} for all.
#' @param gene_id_type Either \code{"symbol"} or \code{"entrez"}: the
#'   identifier used to build the TERM2GENE table and reported in the
#'   \code{EnrichedGenes} column.
#' @param universe Background gene identifiers, or \code{NULL} for all
#'   genes in the gene sets. Converted to match \code{gene_id_type}, so
#'   Entrez IDs may be supplied whichever mode is in use.
#' @param ... Forwarded to \code{\link{ora}} (\code{minGSSize},
#'   \code{maxGSSize}).
#'
#' @return A data frame of ORA enrichment results.
#' @keywords internal
.run_ora <- function(x, chemicals_meta, cache_dir, pAdjustMethod,
    interaction_types = NULL, gene_id_type = "symbol", universe = NULL, ...) {
    if (gene_id_type == "entrez") {
        if (!is.null(interaction_types)) {
            entrez_list <- .filter_gene_sets(cache_dir, interaction_types)$entrez
        } else {
            entrez_list <- .ctd_cache_load(
                .ctd_bfc(cache_dir), "ChemicalName_GeneEntrezIds")
        }
        term2gene   <- do.call(rbind, lapply(names(entrez_list), function(chem)
            data.frame(term = chem, gene = entrez_list[[chem]],
                       stringsAsFactors = FALSE)))
        input_genes <- as.character(x$EntrezID)
    } else {
        if (!is.null(interaction_types)) {
            gs       <- .filter_gene_sets(cache_dir, interaction_types)
            term2gene <- gs$symbols
        } else {
            term2gene <- .ctd_cache_load(
                .ctd_bfc(cache_dir), "ChemicalName_GeneSymbols")
        }
        sym_map <- suppressMessages(AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys = as.character(x$EntrezID),
            column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
        ))
        # Fallback to Entrez ID for unmapped genes
        input_genes <- ifelse(is.na(sym_map), names(sym_map), sym_map)
        # The universe has to go through the same conversion as the input.
        # Left as Entrez IDs it would not intersect a symbol-keyed
        # background at all, and an empty background silently produces an
        # empty result rather than an error.
        universe <- .to_symbols(universe)
    }

    res <- ora(
        term2gene, input_genes,
        pAdjustMethod = pAdjustMethod,
        universe = universe,
        ...
    )

    .format_enrichment_result(res, chemicals_meta, pAdjustMethod,
        method = "ORA",
        rename = c(
            pvalue         = "PValue",
            BgRatio        = "BackgroundRatio",
            geneID         = "EnrichedGenes",
            foldEnrichment = "FoldEnrichment"
        ),
        drop = c("p.adjust")
    )
}

#' Convert gene identifiers to HGNC symbols where possible
#'
#' Entrez IDs that do not map keep their original value, matching how the
#' input gene list is handled, so a partially mappable universe narrows
#' the background rather than emptying it. Values that are already
#' symbols pass through: they simply do not map, and are kept.
#'
#' @param ids Vector of gene identifiers, or \code{NULL}.
#' @return A character vector of symbols, or \code{NULL} for \code{NULL}
#'   input.
#' @keywords internal
.to_symbols <- function(ids) {
    if (is.null(ids)) return(NULL)
    ids <- as.character(ids)
    # Only all-digit values can be Entrez IDs. Asking AnnotationDbi to map
    # a vector of symbols is not merely useless, it errors when none of
    # the keys is valid, so a universe already given as symbols would
    # bring the call down.
    is_entrez <- grepl("^[0-9]+$", ids)
    if (!any(is_entrez)) return(ids)
    mapped <- suppressMessages(suppressWarnings(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db, keys = ids[is_entrez],
        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
    )))
    ids[is_entrez] <- ifelse(is.na(mapped), ids[is_entrez], mapped)
    unname(ids)
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
#' @param interaction_types Character vector of CTD \code{InteractionActions}
#'   values to retain when building gene sets, or \code{NULL} for all.
#' @param gene_id_type Either \code{"symbol"} or \code{"entrez"}: the
#'   identifier reported in the \code{EnrichedGenes} column.
#' @param ... Forwarded to \code{\link{gsea}} (e.g. \code{minSize},
#'   \code{maxSize}).
#'
#' @return A data frame of GSEA enrichment results, formatted by
#'   \code{\link{.format_enrichment_result}}.
#' @keywords internal
.run_gsea <- function(x, chemicals_meta, cache_dir, pAdjustMethod,
    interaction_types = NULL, gene_id_type = "symbol", ...) {
    if (!is.null(interaction_types)) {
        entrez_sets <- .filter_gene_sets(cache_dir, interaction_types)$entrez
    } else {
        entrez_sets <- .ctd_cache_load(
            .ctd_bfc(cache_dir), "ChemicalName_GeneEntrezIds")
    }
    gene_table <- as.data.frame(x)
    gene_table <- gene_table[!is.na(gene_table$EntrezID), ]

    # Build the label column used by .annotate_genes() for EnrichedGenes
    if (gene_id_type == "entrez") {
        gene_table$GeneLabel <- as.character(gene_table$EntrezID)
    } else {
        sym_map <- suppressMessages(AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys = as.character(gene_table$EntrezID),
            column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
        ))
        gene_table$GeneLabel <- ifelse(is.na(sym_map),
                                       as.character(gene_table$EntrezID),
                                       sym_map)
    }
    res <- gsea(entrez_sets, gene_table, ...)

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

#' Filter cached gene sets by interaction type
#'
#' Loads \code{ctd_interactions.rda} and rebuilds Entrez-ID gene set lists
#' and symbol TERM2GENE tables retaining only interactions whose
#' \code{InteractionActions} field matches at least one of the requested types.
#'
#' @param cache_dir Directory holding cached CTD \code{.rda} files.
#' @param interaction_types Character vector of interaction types to retain.
#' @return A list with \code{entrez} (named list: ChemicalID → Entrez IDs)
#'   and \code{symbols} (data frame with columns \code{term}, \code{gene}).
#' @keywords internal
.filter_gene_sets <- function(cache_dir, interaction_types) {
    bfc <- .ctd_bfc(cache_dir)
    if (!.ctd_cache_has(bfc, "ctd_interactions"))
        stop(
            "ctd_interactions not found in cache. ",
            "Please re-run import_CTD() to rebuild the cache with ",
            "interaction-type support.",
            call. = FALSE
        )
    ia   <- .ctd_cache_load(bfc, "ctd_interactions")
    keep <- !is.na(ia$InteractionActions) &
        vapply(ia$InteractionActions, function(x)
            any(interaction_types %in% strsplit(x, "|", fixed = TRUE)[[1]]),
            logical(1))
    ia <- ia[keep, ]
    if (nrow(ia) == 0L)
        stop("No interactions remain after filtering by interaction_types. ",
             "Check that the values match the CTD vocabulary.",
             call. = FALSE)
    message("interaction_types filter: retained ", nrow(ia),
            " (ChemicalID, gene) pairs matching [",
            paste(interaction_types, collapse = ", "), "]")

    entrez <- split(ia$EntrezID, ia$ChemicalID)
    entrez <- lapply(entrez, unique)

    all_entrez <- unique(unlist(entrez, use.names = FALSE))
    sym_vec <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = all_entrez, column = "SYMBOL",
        keytype = "ENTREZID", multiVals = "first"
    ))
    symbols <- do.call(rbind, lapply(names(entrez), function(chem) {
        syms <- sym_vec[entrez[[chem]]]
        syms <- syms[!is.na(syms)]
        if (!length(syms)) return(NULL)
        data.frame(term = chem, gene = unname(syms),
                   stringsAsFactors = FALSE)
    }))
    list(entrez = entrez, symbols = symbols)
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
