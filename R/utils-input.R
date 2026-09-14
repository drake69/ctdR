#' @title Internal input helpers for ctdR
#'
#' @description
#' Shared helpers used by the matrix-based enrichment methods
#' (CAMERA, GSVA): identifier type detection from rownames and gene-set
#' loading from the CTD cache.
#'
#' @keywords internal
#' @name utils-input
NULL

#' Detect identifier type from a character vector
#'
#' @param ids Character vector of identifiers (typically \code{rownames(x)}).
#' @return Either \code{"entrez"} (if every non-empty id is purely numeric)
#'   or \code{"symbol"} otherwise.
#' @keywords internal
.detect_id_type <- function(ids) {
    ids <- as.character(ids)
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) == 0) {
        stop("Cannot auto-detect 'id_type': rownames(x) are empty or all NA.",
            call. = FALSE)
    }
    if (all(grepl("^[0-9]+$", ids))) "entrez" else "symbol"
}

#' Load CTD gene-set list from the user cache
#'
#' @param id_type Either \code{"entrez"} or \code{"symbol"}.
#' @param cache_dir Directory holding the cached \code{.rda} files.
#' @return A named list where names are CTD chemical IDs and elements are
#'   character vectors of gene IDs (Entrez or HGNC SYMBOL).
#' @keywords internal
.load_geneset_list <- function(id_type, cache_dir) {
    id_type <- match.arg(id_type, c("entrez", "symbol"))
    bfc <- .ctd_bfc(cache_dir)
    if (id_type == "entrez") {
        gs <- .ctd_cache_load(bfc, "ChemicalName_GeneEntrezIds")
        gs <- lapply(gs, as.character)
    } else {
        df <- .ctd_cache_load(bfc, "ChemicalName_GeneSymbols")
        df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
        gs <- split(as.character(df$gene), as.character(df$term))
    }
    gs
}

#' Validate that an expression matrix has usable rownames
#'
#' @param x A numeric matrix.
#' @keywords internal
.validate_expr_matrix <- function(x) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("'x' must be a numeric matrix (genes x samples).",
            call. = FALSE)
    }
    .validate_expr_rownames(rownames(x))
    invisible(TRUE)
}

#' Test whether an object is a SummarizedExperiment
#'
#' Uses \code{methods::is()} so the check works without attaching
#' \pkg{SummarizedExperiment}, and keeps the dependency confined to the
#' matrix-based branches.
#'
#' @param x Any object.
#' @return \code{TRUE} if \code{x} inherits from
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @keywords internal
.is_se <- function(x) {
    methods::is(x, "SummarizedExperiment")
}

#' Resolve an assay selector to a positional index
#'
#' @param assays_list The \code{assays()} list of a
#'   \code{SummarizedExperiment}.
#' @param assay A single assay name, a positive integer index, or
#'   \code{NULL} for the first assay.
#' @return An integer index into \code{assays_list}.
#' @keywords internal
.assay_index <- function(assays_list, assay) {
    n <- length(assays_list)
    if (n == 0L) {
        stop("'x' is a SummarizedExperiment with no assays.", call. = FALSE)
    }
    if (is.null(assay)) {
        return(1L)
    }
    if (length(assay) != 1L || is.na(assay)) {
        stop("'assay' must be a single assay name or positive index.",
            call. = FALSE)
    }
    if (is.numeric(assay)) {
        if (assay < 1 || assay > n) {
            stop("'assay' index ", assay, " is out of range: 'x' has ",
                n, " assay(s).", call. = FALSE)
        }
        return(as.integer(assay))
    }
    idx <- match(as.character(assay), names(assays_list))
    if (is.na(idx)) {
        available <- if (is.null(names(assays_list))) {
            "none (assays are unnamed)"
        } else {
            paste(names(assays_list), collapse = ", ")
        }
        stop("assay '", assay, "' not found in 'x'. Available: ", available,
            call. = FALSE)
    }
    idx
}

#' Validate gene identifiers taken from expression rownames
#'
#' @param rn The \code{rownames()} of a matrix or SummarizedExperiment.
#' @keywords internal
.validate_expr_rownames <- function(rn) {
    if (is.null(rn)) {
        stop("'x' must have rownames (gene identifiers).", call. = FALSE)
    }
    if (anyDuplicated(rn)) {
        stop("'x' has duplicated rownames; gene identifiers must be unique.",
            call. = FALSE)
    }
    invisible(TRUE)
}

#' Select one assay of a SummarizedExperiment, keeping the container
#'
#' Reduces \code{x} to the single requested assay, which then becomes the
#' first (and only) one. Methods that consume a
#' \code{SummarizedExperiment} directly, such as
#' \code{\link[GSVA]{gsvaParam}}, pick the first assay by default, so this
#' avoids having to pass an assay name downstream and works whether or not
#' the assays are named.
#'
#' @param x A \code{SummarizedExperiment}.
#' @param assay Assay name, index, or \code{NULL} for the first.
#' @return \code{x} with a single assay, validated.
#' @importFrom SummarizedExperiment assays
#' @keywords internal
.select_se_assay <- function(x, assay = NULL) {
    a <- SummarizedExperiment::assays(x)
    idx <- .assay_index(a, assay)
    SummarizedExperiment::assays(x) <- a[idx]
    .validate_expr_rownames(rownames(x))
    x
}

#' Coerce an expression input to a plain numeric matrix
#'
#' Single coercion point for the matrix-based methods. Accepts either a
#' numeric matrix (returned unchanged) or a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, from which the
#' requested assay is extracted.
#'
#' @param x A numeric matrix (genes x samples) or a
#'   \code{SummarizedExperiment}.
#' @param assay Assay name, index, or \code{NULL} for the first. Ignored
#'   when \code{x} is already a matrix.
#' @return A validated numeric matrix with genes in rows.
#' @keywords internal
.as_expr_matrix <- function(x, assay = NULL) {
    if (.is_se(x)) {
        a <- SummarizedExperiment::assays(x)
        idx <- .assay_index(a, assay)
        x <- as.matrix(a[[idx]])
    }
    .validate_expr_matrix(x)
    x
}
