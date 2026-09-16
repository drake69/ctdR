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
#'   \item Saves four cached objects via \pkg{BiocFileCache}:
#'     \code{chemicals}, \code{ChemicalName_GeneEntrezIds},
#'     \code{ChemicalName_GeneSymbols}, and \code{ctd_interactions}
#'     (a long-format table of chemical--gene--action triples used by
#'     \code{enrichment_CTD(interaction_types = ...)}).
#' }
#'
#' The cache is managed by \pkg{BiocFileCache} under
#' \code{tools::R_user_dir("ctdR", "cache")}. To re-import (e.g. after
#' downloading a newer CTD release), simply call \code{import_CTD()} again —
#' existing cache resources are overwritten.
#'
#' Filtering by interaction type (e.g., to retain only
#' \code{"increases^expression"} interactions) is done at enrichment time via
#' the \code{interaction_types} argument of \code{\link{enrichment_CTD}},
#' not here — so a single import supports any combination of filters without
#' re-running this step.
#'
#' @param file_path Character. A local path \strong{or} a URL to the CTD
#'   chemical-gene interactions file (\code{CTD_chem_gene_ixns.csv} or
#'   \code{CTD_chem_gene_ixns.csv.gz}). Remote URLs are downloaded and cached
#'   with \pkg{BiocFileCache}; the package assumes no default URL, so you supply
#'   the source and a one-time CTD data-licensing reminder is shown on first
#'   remote fetch.
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
#' @importFrom BiocFileCache bfccache
#' @export
import_CTD <- function(file_path) {
    src <- .resolve_ctd_source(file_path)

    bfc <- .ctd_bfc()

    t0 <- proc.time()[["elapsed"]]

    CTD_chem_gene_ixns <- .read_and_validate_ctd(src)

    chemicals_ids <- unique(CTD_chem_gene_ixns$ChemicalID)
    gene_maps <- .map_chemical_genes(CTD_chem_gene_ixns, chemicals_ids)
    chemicals  <- .deduplicate_chemicals(
        CTD_chem_gene_ixns[, c("ChemicalID", "ChemicalName")]
    )

    interactions <- .build_interaction_table(CTD_chem_gene_ixns)

    hdr <- attr(CTD_chem_gene_ixns, "ctd_header")
    provenance <- .ctd_provenance_new(
        report_created = hdr$report_created,
        # The source as the user gave it, not the resolved path: for a URL
        # the resolved path is a cache filename that names nothing.
        source = file_path,
        n_chemicals = nrow(chemicals),
        n_interactions = nrow(interactions)
    )

    .save_ctd_cache(bfc, chemicals,
                    gene_maps$entrez, gene_maps$symbols, interactions,
                    provenance)

    elapsed <- proc.time()[["elapsed"]] - t0
    release <- if (is.na(provenance$report_created))
        "not stated in the file" else provenance$report_created
    message(sprintf(
        "CTD data cached successfully in: %s", BiocFileCache::bfccache(bfc)))
    message(sprintf("  %d chemicals | %d unique genes | %.0f s",
        nrow(chemicals),
        length(unique(interactions$EntrezID)),
        elapsed))
    message("  CTD release: ", release)
    invisible(NULL)
}

#' Deduplicate the chemical metadata table
#'
#' Ensures one row per \code{ChemicalID}. Warns if the same ID appears with
#' multiple names (CTD data quality issue — first occurrence is kept). Reports
#' as a message if the same name is shared by multiple IDs (legitimate:
#' parent compound and its derivatives each have distinct MeSH IDs).
#'
#' @param chem_df Data frame with columns \code{ChemicalID} and
#'   \code{ChemicalName}.
#' @return A data frame with one row per unique \code{ChemicalID}.
#' @keywords internal
.deduplicate_chemicals <- function(chem_df) {
    dup_id <- duplicated(chem_df$ChemicalID)
    if (any(dup_id)) {
        ambiguous <- unique(chem_df$ChemicalID[dup_id])
        warning(
            length(ambiguous), " ChemicalID(s) appear with more than one ",
            "ChemicalName in the CTD file; only the first name per ID is ",
            "retained. Affected IDs: ",
            paste(ambiguous[seq_len(min(5L, length(ambiguous)))],
                collapse = ", "),
            if (length(ambiguous) > 5L)
                paste0(" ... (and ", length(ambiguous) - 5L, " more)"),
            call. = FALSE
        )
    }
    out <- chem_df[!dup_id, ]

    dup_name <- unique(out$ChemicalName[duplicated(out$ChemicalName)])
    if (length(dup_name) > 0L)
        message(
            length(dup_name), " ChemicalName(s) are shared by more than one ",
            "ChemicalID (e.g., a parent compound and its derivatives). ",
            "Each ChemicalID is treated as a distinct gene set."
        )

    rownames(out) <- NULL
    out
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

#' Column names a CTD chemical-gene interactions file must provide
#'
#' Used both to locate the column-name line inside the file header and to
#' validate the result of reading it.
#' @return Character vector of expected CTD field names.
#' @keywords internal
.ctd_known_fields <- function() {
    c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID", "Interaction",
        "InteractionActions", "PubMedIDs")
}

#' Parse the comment header of a CTD file
#'
#' A CTD download has no header row. The field names sit \emph{inside} the
#' commented preamble, on the line after \code{# Fields:}, and the release
#' date sits on the line beginning \code{# Report created:}.
#'
#' The column-name line is located by looking for the commented line that
#' lists at least \code{min_match} of the names in
#' \code{\link{.ctd_known_fields}}, rather than by counting header lines.
#' Counting is the more fragile assumption of the two: the package already
#' depends on those names everywhere (\code{OrganismID}, \code{ChemicalID}
#' and the rest are referenced throughout), so matching them adds no new
#' dependency, whereas a hard-coded line count adds one that buys nothing.
#' The failure modes differ too. If CTD inserts a comment line, a fixed
#' count shifts every column silently and the analysis proceeds on
#' misaligned data; if CTD renames a column, this search finds nothing and
#' stops before a single record is read.
#'
#' @param file_path Path to the CTD CSV file.
#' @param min_match Integer. How many known field names a commented line
#'   must list to be taken as the column-name line. Three, so that a
#'   passing mention of one name in the licence text cannot be mistaken
#'   for the real thing.
#' @param n_peek Integer. How many lines to read from the top of the file.
#'   The header of a CTD download is 29 lines, so 50 leaves room for it to
#'   grow while keeping this a single bounded read of a file that is
#'   hundreds of megabytes on disk.
#'
#' @return A list with \code{fields} (character vector of column names in
#'   file order), \code{report_created} (the CTD release string, or
#'   \code{NA_character_} when the file does not carry one) and
#'   \code{n_comment_lines}.
#' @keywords internal
.parse_ctd_header <- function(file_path, min_match = 3L, n_peek = 50L) {
    known <- .ctd_known_fields()
    peek <- readLines(file_path, n = n_peek, warn = FALSE)
    if (!length(peek))
        stop("'", file_path, "' is empty.", call. = FALSE)
    # Only the leading run of comment lines is the header; anything after
    # the first data row is data and must not be searched for field names.
    is_comment <- grepl("^#", peek)
    first_data <- match(FALSE, is_comment, nomatch = length(peek) + 1L)
    hdr <- peek[seq_len(first_data - 1L)]

    split_fields <- function(line)
        trimws(strsplit(sub("^#[[:space:]]*", "", line), ",",
            fixed = TRUE)[[1]])
    n_known <- vapply(hdr,
        function(l) sum(split_fields(l) %in% known), integer(1),
        USE.NAMES = FALSE)

    # Among the lines that match best, take the LAST. The field-name line
    # sits immediately before the data by convention, so anything earlier
    # naming the same columns is prose in the preamble, or a duplicate.
    best_i <- if (length(n_known))
        max(which(n_known == max(n_known))) else NA_integer_

    if (!length(n_known) || max(n_known) < min_match) {
        best <- if (length(n_known)) hdr[best_i] else NA_character_
        stop("Could not find the column-name line in the header of '",
            file_path, "'.\n",
            "  Expected: a commented line listing at least ", min_match,
            " of these names, separated by commas:\n    ",
            paste(known, collapse = ", "), "\n",
            "  Scanned ", length(hdr), " commented line(s)",
            if (length(hdr) >= n_peek)
                paste0(" (the first ", n_peek, " lines were all comments; ",
                    "raise 'n_peek' if this file has a longer header)") else "",
            "; the closest listed ",
            if (length(n_known)) max(n_known) else 0L, " of them",
            if (!is.na(best)) paste0(":\n    ", substr(best, 1, 120)) else ".",
            call. = FALSE
        )
    }

    created <- grep("^#[[:space:]]*Report created:", hdr, value = TRUE)
    list(
        fields = split_fields(hdr[[best_i]]),
        report_created = if (length(created))
            trimws(sub("^#[[:space:]]*Report created:[[:space:]]*", "",
                created[[1]])) else NA_character_,
        n_comment_lines = length(hdr)
    )
}

#' Read and validate a CTD CSV file
#' @param file_path Path to the CTD CSV file.
#' @return A filtered data frame of human CTD interactions, carrying the
#'   parsed header in its \code{"ctd_header"} attribute.
#' @keywords internal
.read_and_validate_ctd <- function(file_path) {
    message("Reading CTD chemical-gene interactions from: ", file_path)
    hdr <- .parse_ctd_header(file_path)
    # comment = "#" rather than a hard-coded skip, as requested in review.
    # Note for future maintenance: readr drops everything after a "#"
    # anywhere in a line, not only at the start. No record in the CTD
    # chemical-gene file carries one, so this is safe here; it would stop
    # being safe if CTD ever admitted "#" into a field.
    ctd <- suppressWarnings(readr::read_csv(file_path,
        comment = "#", col_names = hdr$fields, show_col_types = FALSE,
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
    ctd <- subset(ctd, ctd$OrganismID == 9606)
    message("Filtered to ", nrow(ctd), " human interactions")
    attr(ctd, "ctd_header") <- hdr
    ctd
}

#' Save CTD cache files
#' @param bfc A \code{BiocFileCache} object (from \code{.ctd_bfc()}).
#' @param chemicals Data frame of chemical IDs and names.
#' @param entrez Named list of Entrez IDs per chemical.
#' @param symbols Data frame of term-gene symbol mappings.
#' @param interactions Long-format data frame (ChemicalID, EntrezID,
#'   InteractionActions).
#' @param provenance A \code{ctd_provenance} record, or \code{NULL} to skip
#'   writing one (which is what a caller that has none should pass).
#' @return Invisible \code{NULL}.
#' @keywords internal
.save_ctd_cache <- function(bfc, chemicals, entrez, symbols, interactions,
    provenance = NULL) {
    .ctd_cache_save(bfc, "chemicals", chemicals)
    .ctd_cache_save(bfc, "ChemicalName_GeneEntrezIds", entrez)
    .ctd_cache_save(bfc, "ChemicalName_GeneSymbols", symbols)
    .ctd_cache_save(bfc, "ctd_interactions", interactions)
    if (!is.null(provenance))
        .ctd_cache_save(bfc, "ctd_provenance", provenance)
    invisible(NULL)
}
