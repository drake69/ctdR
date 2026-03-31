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
#'   \item Saves three cached objects (\code{chemicals}, \code{ChemicalName_GeneEntrezIds},
#'         \code{ChemicalName_GeneSymbols}) to the user cache directory
#'         (see \code{\link[rappdirs]{user_cache_dir}}).
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
#' # Check that import_CTD validates the file path:
#' file.exists("/nonexistent/path.csv")  # FALSE
#'
#' \donttest{
#' # Full workflow (requires downloaded CTD data):
#' # 1. Download CTD_chem_gene_ixns.csv.gz from:
#' #    https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz
#' # 2. Decompress: gunzip CTD_chem_gene_ixns.csv.gz
#' # 3. Import:
#' import_CTD("~/Downloads/CTD_chem_gene_ixns.csv")
#' }
#'
#' @importFrom rappdirs user_cache_dir
#' @export
import_CTD <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path, "\n",
         "Please download CTD_chem_gene_ixns.csv.gz from:\n",
         "  https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz\n",
         "Decompress it (gunzip CTD_chem_gene_ixns.csv.gz) and provide the path to the .csv file.\n",
         call. = FALSE)
  }

  cache_dir <- rappdirs::user_cache_dir("ctdR")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  message("Reading CTD chemical-gene interactions from: ", file_path)
  CTD_chem_gene_ixns <- readr::read_csv(file_path, skip = 27,
                                         show_col_types = FALSE)

  if (nrow(CTD_chem_gene_ixns) < 2) {
    stop("The file appears to be empty or does not contain enough data rows. ",
         "Ensure you are using the correct CTD file: CTD_chem_gene_ixns.csv",
         call. = FALSE)
  }

  # Validate expected CTD columns are present
  required_cols <- c("ChemicalID", "CasRN", "GeneSymbol", "GeneID",
                     "GeneForms", "Organism", "OrganismID", "Interaction",
                     "InteractionActions", "PubMedIDs")
  found_cols <- colnames(CTD_chem_gene_ixns)
  missing_cols <- setdiff(required_cols, found_cols)
  if (length(missing_cols) > 0) {
    stop("The input file does not match the expected CTD chemical-gene ",
         "interactions format.\n",
         "  Missing columns: ", paste(missing_cols, collapse = ", "), "\n",
         "  Found columns: ", paste(found_cols, collapse = ", "), "\n",
         "Please ensure you downloaded 'CTD_chem_gene_ixns.csv.gz' from:\n",
         "  https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz",
         call. = FALSE)
  }

  CTD_chem_gene_ixns <- CTD_chem_gene_ixns[-1, ]

  # rename first column to ChemicalName
  names(CTD_chem_gene_ixns)[1] <- "ChemicalName"

  # Filter only human (OrganismID 9606)
  CTD_chem_gene_ixns <- subset(CTD_chem_gene_ixns, CTD_chem_gene_ixns$OrganismID == 9606)
  message("Filtered to ", nrow(CTD_chem_gene_ixns), " human interactions")

  chemicals_ids <- unique(CTD_chem_gene_ixns$ChemicalID)

  ChemicalName_GeneEntrezIds <- list()
  ChemicalName_GeneSymbols <- data.frame()

  message("Mapping genes for ", length(chemicals_ids), " chemicals (this may take a while)...")
  for (i in seq_along(chemicals_ids)) {
    chemical <- chemicals_ids[i]
    gene_entrezids <- unique(CTD_chem_gene_ixns[CTD_chem_gene_ixns$ChemicalID == chemical, "GeneID"])
    gene_entrezids <- as.data.frame(gene_entrezids)[, 1]
    ChemicalName_GeneEntrezIds[[chemical]] <- gene_entrezids

    gene_entrezids <- as.character(gene_entrezids)
    tryCatch({
      gene_symbols <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = gene_entrezids,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first")
      tt <- expand.grid(term = chemical, gene = gene_symbols)
      ChemicalName_GeneSymbols <- plyr::rbind.fill(tt, ChemicalName_GeneSymbols)
    }, error = function(e) {
      # skip chemicals with unmappable genes
    })
  }

  chemicals <- unique(CTD_chem_gene_ixns[, c("ChemicalID", "ChemicalName")])

  save(chemicals, file = file.path(cache_dir, "chemicals.rda"))
  save(ChemicalName_GeneEntrezIds, file = file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"))
  save(ChemicalName_GeneSymbols, file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"))

  message("CTD data cached successfully in: ", cache_dir)
  invisible(NULL)
}
