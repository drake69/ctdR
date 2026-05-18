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
    if (id_type == "entrez") {
        e <- new.env(parent = emptyenv())
        load(file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"),
            envir = e)
        gs <- e$ChemicalName_GeneEntrezIds
        gs <- lapply(gs, as.character)
    } else {
        e <- new.env(parent = emptyenv())
        load(file.path(cache_dir, "ChemicalName_GeneSymbols.rda"),
            envir = e)
        df <- e$ChemicalName_GeneSymbols
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
    if (is.null(rownames(x))) {
        stop("'x' must have rownames (gene identifiers).", call. = FALSE)
    }
    if (anyDuplicated(rownames(x))) {
        stop("'x' has duplicated rownames; gene identifiers must be unique.",
            call. = FALSE)
    }
    invisible(TRUE)
}
