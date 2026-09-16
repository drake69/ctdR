#' @title CTD Chemical-Gene Interactions File
#'
#' @description
#' \code{CTDFile} is a \code{\link[BiocIO]{BiocFile}} subclass representing a
#' Comparative Toxicogenomics Database (CTD) chemical-gene interactions file
#' (\code{CTD_chem_gene_ixns.csv} or its gzip-compressed \code{.csv.gz} form).
#' Pair it with \code{\link[BiocIO]{import}} to read the file into R following
#' the Bioconductor \pkg{BiocIO} import/export convention, as an alternative to
#' the side-effecting \code{\link{import_CTD}} caching workflow.
#'
#' In \pkg{ctdR}, \strong{CTD} always denotes the \emph{Comparative
#' Toxicogenomics Database} (\url{https://ctdbase.org}). The class is named
#' \code{CTDFile} (not \code{CTD}) to avoid confusion with unrelated software
#' that shares the acronym (the CRAN package \pkg{CTD} is a graph algorithm; the
#' \code{ctd} object in \pkg{EWCE} is a cell-type dataset).
#'
#' @section Data Licensing Disclaimer:
#' This package does \strong{not} bundle or redistribute any CTD data. The
#' Comparative Toxicogenomics Database is maintained by NC State University and
#' its data are subject to specific licensing terms. Users are responsible for
#' downloading the data directly from \url{https://ctdbase.org} and for
#' complying with the CTD Terms of Service
#' (\url{https://ctdbase.org/about/legal.jsp}).
#'
#' @param resource Character scalar: a local path or a URL to the CTD
#'   chemical-gene interactions file.
#' @param con A \code{CTDFile} object to import.
#' @param format Ignored; present for compatibility with the \pkg{BiocIO}
#'   \code{import} generic.
#' @param text Ignored; present for compatibility with the \pkg{BiocIO}
#'   \code{import} generic.
#' @param ... Additional arguments (currently unused).
#'
#' @return \code{CTDFile()} returns a \code{CTDFile} object. \code{import()} on a
#'   \code{CTDFile} returns a \code{\link[S4Vectors]{DataFrame}} of validated
#'   human (OrganismID 9606) CTD chemical-gene interactions.
#'
#' @seealso \code{\link{import_CTD}} for the caching workflow required before
#'   \code{\link{enrichment_CTD}}; \code{\link[BiocIO]{BiocFile}} for the parent
#'   class.
#'
#' @examples
#' sample_file <- system.file(
#'     "extdata", "CTD_chem_gene_ixns_sample.csv",
#'     package = "ctdR"
#' )
#' ctd_file <- CTDFile(sample_file)
#' ctd_file
#' interactions <- BiocIO::import(ctd_file)
#' head(interactions)
#'
#' @name CTDFile
#' @rdname CTDFile
#' @aliases CTDFile-class
#' @importClassesFrom BiocIO BiocFile
#' @importMethodsFrom BiocIO import
#' @importFrom BiocIO resource
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#' @exportClass CTDFile
.CTDFile <- setClass("CTDFile", contains = "BiocFile")

#' @rdname CTDFile
#' @export
CTDFile <- function(resource) {
    if (!is.character(resource) || length(resource) != 1L || is.na(resource)) {
        stop("'resource' must be a single non-NA character path or URL.",
            call. = FALSE
        )
    }
    .CTDFile(resource = resource)
}

#' @rdname CTDFile
#' @exportMethod import
setMethod("import", c(con = "CTDFile"), function(con, format, text, ...) {
    src <- .resolve_ctd_source(resource(con))
    ctd <- .read_and_validate_ctd(src)
    hdr <- attr(ctd, "ctd_header")
    out <- S4Vectors::DataFrame(
        as.data.frame(ctd, stringsAsFactors = FALSE),
        check.names = FALSE
    )
    # A DataFrame has a metadata() slot, so the CTD release travels with
    # the object read here just as it does with enrichment results.
    .attach_provenance(out, .ctd_provenance_new(
        report_created = hdr$report_created,
        source = resource(con),
        n_chemicals = length(unique(ctd$ChemicalID)),
        n_interactions = nrow(ctd)
    ))
})
