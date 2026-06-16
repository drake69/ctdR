#' @title Import CTD Chemical-Gene Interaction Data
#'
#' @description
#' Parses, filters, and caches the CTD chemical-gene interactions file so that
#' it can be used by \code{\link{enrichment_CTD}}. This function must be called
#' \strong{once} before running any enrichment analysis.
#'
#' The raw data file must be downloaded manually from the CTD website. The
#' required file is \strong{CTD_chem_gene_ixns.csv.gz}, available at
#' \url{https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz}.
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
#'   \item Saves four cached objects to the user cache directory:
#'     \code{chemicals}, \code{ChemicalName_GeneEntrezIds},
#'     \code{ChemicalName_GeneSymbols}, and \code{ctd_interactions}
#'     (a long-format table of chemical--gene--action triples used by
#'     \code{enrichment_CTD(interaction_types = ...)}).
#' }
#'
#' The cache is stored under \code{rappdirs::user_cache_dir("ctdR")}. To
#' re-import (e.g. after downloading a newer CTD release), simply call
#' \code{import_CTD()} again — existing cache files will be overwritten.
#'
#' Filtering by interaction type (e.g., to retain only
#' \code{"increases^expression"} interactions) is done at enrichment time via
#' the \code{interaction_types} argument of \code{\link{enrichment_CTD}},
#' not here — so a single import supports any combination of filters without
#' re-running this step.
#'
#' @param file_path Character. Path to the CTD chemical-gene interactions file
#'   (\code{CTD_chem_gene_ixns.csv} or \code{CTD_chem_gene_ixns.csv.gz}).
#'
#' @return Invisible \code{NULL}. Called for its side effect of caching the
#'   processed data.
#'
#' @seealso \code{\link{enrichment_CTD}} for running enrichment analysis after
#'   import, including the \code{interaction_types} filter.
#'
#' @examples
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
            call. = FALSE
        )
    }

    cache_dir <- rappdirs::user_cache_dir("ctdR")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    CTD_chem_gene_ixns <- .read_and_validate_ctd(file_path)

    chemicals_ids <- unique(CTD_chem_gene_ixns$ChemicalID)
    gene_maps <- .map_chemical_genes(CTD_chem_gene_ixns, chemicals_ids)
    chemicals <- unique(CTD_chem_gene_ixns[, c("ChemicalID", "ChemicalName")])

    interactions <- .build_interaction_table(CTD_chem_gene_ixns)

    .save_ctd_cache(cache_dir, chemicals,
                    gene_maps$entrez, gene_maps$symbols, interactions)
    message("CTD data cached successfully in: ", cache_dir)
    invisible(NULL)
}

#' Build the long-format chemical–gene–action interaction table
#'
#' For each (ChemicalID, GeneID) pair, collapses all observed
#' InteractionActions values (pipe-separated) into a single string.
#' This table is cached and used at enrichment time when
#' \code{interaction_types} is specified.
#'
#' @param ctd_data Filtered CTD interaction data frame.
#' @return A data frame with columns \code{ChemicalID}, \code{EntrezID},
#'   and \code{InteractionActions} (pipe-collapsed per pair).
#' @keywords internal
.build_interaction_table <- function(ctd_data) {
    pairs <- ctd_data[, c("ChemicalID", "GeneID", "InteractionActions")]
    pairs <- pairs[!is.na(pairs$GeneID) & nzchar(pairs$GeneID), ]
    pairs$GeneID <- as.character(pairs$GeneID)
    # Collapse multiple InteractionActions per (ChemicalID, GeneID) pair
    result <- do.call(rbind, lapply(
        split(pairs, paste(pairs$ChemicalID, pairs$GeneID, sep = "\t")),
        function(g) {
            ia_all <- unlist(strsplit(g$InteractionActions, "|", fixed = TRUE))
            ia_all <- unique(ia_all[!is.na(ia_all) & nzchar(ia_all)])
            data.frame(
                ChemicalID         = g$ChemicalID[1],
                EntrezID           = g$GeneID[1],
                InteractionActions = if (length(ia_all)) paste(ia_all, collapse = "|") else NA_character_,
                stringsAsFactors   = FALSE
            )
        }
    ))
    rownames(result) <- NULL
    result
}

#' Map chemicals to gene IDs and symbols
#' @param ctd_data Filtered CTD interaction data frame.
#' @param chemicals_ids Character vector of chemical IDs.
#' @return A list with \code{entrez} (named list) and
#'   \code{symbols} (data frame).
#' @keywords internal
.map_chemical_genes <- function(ctd_data, chemicals_ids) {
    message("Mapping genes for ", length(chemicals_ids), " chemicals...")

    entrez_map <- lapply(chemicals_ids, function(chem) {
        unique(as.character(
            ctd_data[ctd_data$ChemicalID == chem, "GeneID"][[1]]
        ))
    })
    names(entrez_map) <- chemicals_ids

    all_entrez <- unique(unlist(entrez_map, use.names = FALSE))
    sym_vec <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys      = all_entrez,
        column    = "SYMBOL",
        keytype   = "ENTREZID",
        multiVals = "first"
    ))

    symbols_df <- do.call(rbind, lapply(chemicals_ids, function(chem) {
        syms <- sym_vec[entrez_map[[chem]]]
        syms <- syms[!is.na(syms)]
        if (length(syms) == 0L) return(NULL)
        data.frame(term = chem, gene = unname(syms),
                   stringsAsFactors = FALSE)
    }))

    list(entrez = entrez_map, symbols = symbols_df)
}

#' Read and validate a CTD CSV file
#' @param file_path Path to the CTD CSV file.
#' @return A filtered data frame of human CTD interactions.
#' @keywords internal
.read_and_validate_ctd <- function(file_path) {
    message("Reading CTD chemical-gene interactions from: ", file_path)
    ctd <- suppressWarnings(readr::read_csv(file_path,
        skip = 27, show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_character())
    ))
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
#' @param interactions Long-format data frame (ChemicalID, EntrezID,
#'   InteractionActions).
#' @return Invisible \code{NULL}.
#' @keywords internal
.save_ctd_cache <- function(cache_dir, chemicals, entrez, symbols,
    interactions) {
    ChemicalName_GeneEntrezIds <- entrez
    ChemicalName_GeneSymbols   <- symbols
    ctd_interactions           <- interactions
    save(chemicals,
        file = file.path(cache_dir, "chemicals.rda"))
    save(ChemicalName_GeneEntrezIds,
        file = file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"))
    save(ChemicalName_GeneSymbols,
        file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"))
    save(ctd_interactions,
        file = file.path(cache_dir, "ctd_interactions.rda"))
}
