#' @title Provenance of the CTD data an analysis ran on
#'
#' @description
#' Returns the record of which CTD release produced a result: the
#' \code{Report created} date CTD stamps into its own file header, where
#' the file came from, when it was imported, and how much of it was kept.
#'
#' An analysis is only reproducible if the version of the data behind it
#' can be named. CTD is re-released continuously and its downloads are not
#' versioned in the filename, so the release date inside the header is the
#' only thing that identifies which snapshot a result came from. This
#' carries it from the file all the way to the object you report on.
#'
#' @details
#' Where the record is stored depends on what the object is, because the
#' methods do not all return the same container. Objects that provide a
#' \code{metadata()} slot keep it there, which covers the
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} returned by
#' GSVA and the \code{DataFrame} returned by importing a
#' \code{\link{CTDFile}}. The data frames returned by ORA, GSEA and
#' CAMERA, and a plain score matrix, have no such slot, so there it rides
#' on an attribute. Use this accessor rather than reaching for either
#' directly, so that code keeps working whichever method produced the
#' object.
#'
#' The record survives subsetting, ordering, \code{head()} and the common
#' \pkg{dplyr} verbs. It does not survive \code{merge()} or
#' \code{subset()}, which drop attributes; retrieve it before those if you
#' need it afterwards.
#'
#' @param x An object returned by \code{\link{enrichment_CTD}}, or the
#'   \code{DataFrame} returned by importing a \code{\link{CTDFile}}.
#'   Omit it to read the record of the data currently cached, which
#'   answers "which release am I about to analyse" before any analysis
#'   has been run.
#'
#' @return An object of class \code{ctd_provenance}: a list with
#'   \code{report_created} (the CTD release string, \code{NA} if the file
#'   carried none), \code{source}, \code{accessed}, \code{n_chemicals},
#'   \code{n_interactions} and \code{ctdR_version}. Returns \code{NULL},
#'   with a warning, when \code{x} carries no record.
#'
#' @seealso \code{\link{import_CTD}}, which reads the record, and
#'   \code{\link{enrichment_CTD}}, which attaches it to its results.
#'
#' @examples
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
#' genes <- data.frame(
#'     EntrezID = c("7124", "3569", "7157", "672", "1956"),
#'     pvalue = c(0.001, 0.003, 0.01, 0.02, 0.05)
#' )
#' # Which release is cached, before running anything:
#' ctd_provenance()
#'
#' res <- enrichment_CTD(genes, method = "ORA")
#' ctd_provenance(res)
#'
#' @export
ctd_provenance <- function(x) {
    if (missing(x)) {
        prov <- .ctd_provenance_cached(.ctd_bfc())
        if (is.null(prov)) {
            warning("No CTD provenance in the cache. Either no data has ",
                "been imported yet, or it was imported by a version of ",
                "ctdR that did not record one: re-run import_CTD().",
                call. = FALSE)
        }
        return(prov)
    }
    prov <- if (.has_metadata_slot(x)) {
        S4Vectors::metadata(x)[["ctd_provenance"]]
    } else {
        attr(x, "ctd_provenance", exact = TRUE)
    }
    if (is.null(prov)) {
        warning("No CTD provenance attached to this object. Operations ",
            "that drop attributes, such as merge() and subset(), lose it; ",
            "retrieve it from the original result.",
            call. = FALSE)
        return(NULL)
    }
    prov
}

#' Build a provenance record
#'
#' @param report_created CTD's own release string, or \code{NA_character_}.
#' @param source The path or URL the user supplied, before resolution.
#' @param n_chemicals Number of chemicals retained after filtering.
#' @param n_interactions Number of chemical-gene pairs retained.
#' @return An object of class \code{ctd_provenance}.
#' @keywords internal
.ctd_provenance_new <- function(report_created, source,
    n_chemicals, n_interactions) {
    structure(
        list(
            report_created = report_created,
            source = source,
            accessed = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            n_chemicals = as.integer(n_chemicals),
            n_interactions = as.integer(n_interactions),
            ctdR_version = as.character(utils::packageVersion("ctdR"))
        ),
        class = "ctd_provenance"
    )
}

#' Attach a provenance record to a result
#'
#' Uses \code{metadata()} on a \code{SummarizedExperiment}, which is the
#' slot the class provides, and an attribute on anything else.
#'
#' @param x The object to annotate.
#' @param prov A \code{ctd_provenance} record, or \code{NULL} to leave
#'   \code{x} untouched.
#' @return \code{x}, annotated.
#' @keywords internal
.attach_provenance <- function(x, prov) {
    if (is.null(prov)) return(x)
    if (.has_metadata_slot(x)) {
        S4Vectors::metadata(x)[["ctd_provenance"]] <- prov
    } else {
        attr(x, "ctd_provenance") <- prov
    }
    x
}

#' Does this object carry a metadata() slot?
#'
#' \code{Annotated} is the S4Vectors virtual class that provides
#' \code{metadata()}. Both \code{SummarizedExperiment} and \code{DataFrame}
#' extend it, so asking about the slot rather than about a specific class
#' keeps the rule to one line: use the slot where there is one, an
#' attribute where there is not.
#'
#' @param x Any object.
#' @return \code{TRUE} when \code{metadata()} applies to \code{x}.
#' @keywords internal
.has_metadata_slot <- function(x) methods::is(x, "Annotated")

#' Read the cached provenance record, if there is one
#'
#' @param bfc A \pkg{BiocFileCache} object.
#' @return A \code{ctd_provenance} record, or \code{NULL} when the cache
#'   predates provenance tracking.
#' @keywords internal
.ctd_provenance_cached <- function(bfc) {
    if (!.ctd_cache_has(bfc, "ctd_provenance")) return(NULL)
    .ctd_cache_load(bfc, "ctd_provenance")
}

#' @param ... Ignored, present for compatibility with the generic.
#' @rdname ctd_provenance
#' @export
print.ctd_provenance <- function(x, ...) {
    cat("CTD provenance\n")
    cat("  Report created: ",
        if (is.na(x$report_created)) "not stated in the file"
        else x$report_created, "\n", sep = "")
    cat("  Source:         ", x$source, "\n", sep = "")
    cat("  Imported:       ", x$accessed, "\n", sep = "")
    cat("  Retained:       ", format(x$n_chemicals, big.mark = ","),
        " chemicals, ", format(x$n_interactions, big.mark = ","),
        " chemical-gene pairs\n", sep = "")
    cat("  ctdR version:   ", x$ctdR_version, "\n", sep = "")
    invisible(x)
}
