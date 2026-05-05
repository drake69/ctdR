#' @title Import CTD Chemical-Gene Interaction Data
#'
#' @description
#' Parses, filters, and caches the CTD chemical-gene interactions file so that
#' it can be used by \code{\link{enrichment_CTD}}. This function must be called
#' \strong{once} before running any enrichment analysis.
#'
#' The raw data file must be downloaded manually from the CTD website. The
#' required file is \strong{CTD_chem_gene_ixns.csv.gz}, available at
#' \url{https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz}. Decompress it
#' (e.g. \code{gunzip CTD_chem_gene_ixns.csv.gz}) and pass the resulting
#' \code{CTD_chem_gene_ixns.csv} file path to this function.
#'
#' @section Data Licensing Disclaimer:
#' This package does \strong{not} bundle or redistribute any CTD data. The
#' Comparative Toxicogenomics Database is maintained by NC State University and
#' its data are subject to specific licensing terms. Users are responsible for
#' downloading the data directly from \url{https://ctdbase.org} and for
#' complying with the CTD Terms of Service
#' (\url{https://ctdbase.org/about/legal.jsp}). By using this function you
#' acknowledge that you have read and accepted those terms.
#'
#' @details
#' Processing steps performed by \code{import_CTD}:
#' \enumerate{
#'   \item Reads the CSV, skipping the 27 CTD header lines.
#'   \item Filters interactions to \strong{Homo sapiens} only (OrganismID 9606).
#'   \item For each chemical, collects the associated Entrez gene IDs.
#'   \item Maps Entrez IDs to HGNC gene symbols via \pkg{org.Hs.eg.db}.
#'   \item Saves three cached objects (\code{chemicals},
#'     \code{ChemicalName_GeneEntrezIds},
#'     \code{ChemicalName_GeneSymbols}) to the user cache
#'     directory (\code{\link[rappdirs]{user_cache_dir}}).
#' }
#'
#' The cache is stored under \code{rappdirs::user_cache_dir("ctdR")}. To
#' re-import (e.g. after downloading a newer CTD release), simply call
#' \code{import_CTD()} again — existing cache files will be overwritten.
#'
#' @param file_path Character. Path to the decompressed
#'   \code{CTD_chem_gene_ixns.csv} file (e.g.
#'   \code{"~/Downloads/CTD_chem_gene_ixns.csv"}).
#'
#' @return Invisible \code{NULL}. Called for its side effect of caching the
#'   processed data.
#'
#' @seealso \code{\link{enrichment_CTD}} for running enrichment analysis after
#'   import.
#'
#' @examples
#' # Import the bundled sample data:
#' sample_file <- system.file(
#'     "extdata", "CTD_chem_gene_ixns_sample.csv",
#'     package = "ctdR"
#' )
#' import_CTD(sample_file)
#'
#' @importFrom rappdirs user_cache_dir
#' @export
import_CTD <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("File not found: ", file_path, "\n",
            "Please download CTD_chem_gene_ixns.csv.gz from:\n",
            "  https://ctdbase.org/reports/",
            "CTD_chem_gene_ixns.csv.gz\n",
            "Decompress and provide the .csv path.\n",
            call. = FALSE
        )
    }

    cache_dir <- rappdirs::user_cache_dir("ctdR")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    CTD_chem_gene_ixns <- .read_and_validate_ctd(file_path)

    chemicals_ids <- unique(CTD_chem_gene_ixns$ChemicalID)
    gene_maps <- .map_chemical_genes(
        CTD_chem_gene_ixns, chemicals_ids
    )
    chemicals <- unique(
        CTD_chem_gene_ixns[, c("ChemicalID", "ChemicalName")]
    )

    .save_ctd_cache(
        cache_dir, chemicals,
        gene_maps$entrez, gene_maps$symbols
    )
    message("CTD data cached successfully in: ", cache_dir)
    invisible(NULL)
}

#' Map chemicals to gene IDs and symbols
#' @param ctd_data Filtered CTD interaction data frame.
#' @param chemicals_ids Character vector of chemical IDs.
#' @return A list with \code{entrez} (named list) and
#'   \code{symbols} (data frame).
#' @keywords internal
.map_chemical_genes <- function(ctd_data, chemicals_ids) {
    entrez_map <- list()
    symbols_df <- data.frame()
    message(
        "Mapping genes for ", length(chemicals_ids),
        " chemicals (this may take a while)..."
    )
    for (i in seq_along(chemicals_ids)) {
        chemical <- chemicals_ids[i]
        rows <- ctd_data$ChemicalID == chemical
        gene_ids <- unique(ctd_data[rows, "GeneID"])
        gene_ids <- as.data.frame(gene_ids)[, 1]
        entrez_map[[chemical]] <- gene_ids
        gene_ids <- as.character(gene_ids)
        tryCatch(
            {
                syms <- AnnotationDbi::mapIds(
                    org.Hs.eg.db::org.Hs.eg.db,
                    keys = gene_ids,
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first"
                )
                tt <- expand.grid(
                    term = chemical, gene = syms,
                    stringsAsFactors = FALSE
                )
                symbols_df <- plyr::rbind.fill(tt, symbols_df)
            },
            error = function(e) {}
        )
    }
    list(entrez = entrez_map, symbols = symbols_df)
}

#' Read and validate a CTD CSV file
#' @param file_path Path to the CTD CSV file.
#' @return A filtered data frame of human CTD interactions.
#' @keywords internal
.read_and_validate_ctd <- function(file_path) {
    message("Reading CTD chemical-gene interactions from: ", file_path)
    ctd <- readr::read_csv(file_path,
        skip = 27, show_col_types = FALSE
    )
    if (nrow(ctd) < 2) {
        stop("File appears empty or has too few rows. ",
            "Use CTD_chem_gene_ixns.csv",
            call. = FALSE
        )
    }
    required_cols <- c(
        "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID",
        "Interaction", "InteractionActions", "PubMedIDs"
    )
    missing <- setdiff(required_cols, colnames(ctd))
    if (length(missing) > 0) {
        stop("Not a valid CTD file. Missing columns: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    ctd <- ctd[-1, ]
    names(ctd)[1] <- "ChemicalName"
    ctd <- subset(ctd, ctd$OrganismID == 9606)
    message("Filtered to ", nrow(ctd), " human interactions")
    ctd
}

#' Save CTD cache files
#' @param cache_dir Path to the cache directory.
#' @param chemicals Data frame of chemical IDs and names.
#' @param entrez Named list of Entrez IDs per chemical.
#' @param symbols Data frame of term-gene symbol mappings.
#' @return Invisible \code{NULL}. Called for its side effect of saving
#'   cache files.
#' @keywords internal
.save_ctd_cache <- function(cache_dir, chemicals,
    entrez, symbols) {
    ChemicalName_GeneEntrezIds <- entrez
    ChemicalName_GeneSymbols <- symbols
    save(chemicals,
        file = file.path(cache_dir, "chemicals.rda")
    )
    save(ChemicalName_GeneEntrezIds,
        file = file.path(
            cache_dir, "ChemicalName_GeneEntrezIds.rda"
        )
    )
    save(ChemicalName_GeneSymbols,
        file = file.path(
            cache_dir, "ChemicalName_GeneSymbols.rda"
        )
    )
}
