#' @title Export CTD Gene Sets for Existing Enrichment Engines
#'
#' @description
#' Returns the cached CTD chemical gene sets as a named list, in the shape
#' expected by third-party Bioconductor enrichment engines --- for example the
#' \code{gs} argument of \code{EnrichmentBrowser::sbea()}, or the input to
#' \code{GSEABase::GeneSetCollection()}. This lets you run the CTD gene sets
#' through existing infrastructure instead of, or alongside,
#' \code{\link{enrichment_CTD}}.
#'
#' @details
#' \pkg{ctdR} builds \emph{chemical-centric} gene sets: each CTD chemical
#' becomes a gene set of its interacting genes. \code{\link{enrichment_CTD}}
#' consumes these internally; this helper exposes the same gene sets so they
#' integrate with the wider ecosystem. Call \code{\link{import_CTD}} once first
#' to populate the cache.
#'
#' The returned representation is deliberately a plain named list (not a bespoke
#' class): it is accepted directly by \code{EnrichmentBrowser::sbea()} and is a
#' one-line conversion away from a \code{GeneSetCollection}. \pkg{ctdR} does not
#' take a dependency on those engines --- you supply them.
#'
#' @param id_type Either \code{"entrez"} (default) or \code{"symbol"}: the
#'   identifier space of the returned genes, to match your expression data or
#'   downstream engine.
#' @param interaction_types Optional character vector of CTD interaction actions
#'   (e.g. \code{"increases^expression"}) to retain. \code{NULL} (default) keeps
#'   all interactions.
#' @param cache_dir Directory holding the cached CTD files. Defaults to the
#'   \pkg{ctdR} user cache populated by \code{\link{import_CTD}}.
#'
#' @return A named list: names are CTD \code{ChemicalID}s, elements are
#'   character vectors of gene identifiers (Entrez IDs or HGNC symbols).
#'
#' @seealso \code{\link{import_CTD}} to populate the cache;
#'   \code{\link{enrichment_CTD}} for the built-in enrichment interface.
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
#' gene_sets <- as_genesets_CTD("entrez")
#' length(gene_sets)
#' gene_sets[[1]]
#'
#' \dontrun{
#' # Hand the CTD gene sets to an existing engine (EnrichmentBrowser):
#' library(EnrichmentBrowser)
#' res <- sbea(method = "camera", se = my_se, gs = as_genesets_CTD("entrez"))
#' gsRanking(res)
#' }
#'
#' @export
as_genesets_CTD <- function(id_type = c("entrez", "symbol"),
    interaction_types = NULL,
    cache_dir = .ctd_cache_dir()) {
    id_type <- match.arg(id_type)

    if (!.ctd_cache_has(.ctd_bfc(cache_dir), "ChemicalName_GeneEntrezIds")) {
        stop("CTD cache not found in '", cache_dir, "'.\n",
            "Run import_CTD() on your CTD_chem_gene_ixns file first.",
            call. = FALSE
        )
    }

    if (is.null(interaction_types)) {
        return(.load_geneset_list(id_type, cache_dir))
    }

    filt <- .filter_gene_sets(cache_dir, interaction_types)
    if (id_type == "entrez") {
        lapply(filt$entrez, as.character)
    } else {
        split(
            as.character(filt$symbols$gene),
            as.character(filt$symbols$term)
        )
    }
}
