# Replicate Yang et al. (2024) with ctdR, then screen the whole CTD catalogue
# with the same genes. Recomputes, on the CTD release closest to the authors',
# the ctdR columns of inst/extdata/yang2024/yang2024_replication.tsv (set size,
# overlap and the authors' p-value on ctdR's sets) and the ORA over every
# chemical that the vignette "Replicating a published CTD analysis" displays.
#
# Yang Y, Zheng T, Tang Q, et al. Developmental dyslexia genes are selectively
# targeted by diverse environmental pollutants. BMC Psychiatry 24, 509 (2024).
# https://doi.org/10.1186/s12888-024-05952-4 (CC BY 4.0)
#
# CTD distributes only its current release. The release closest to the one
# the authors used ("2023 data") survives as an Internet Archive capture; the
# script downloads it and checks it against the archive's own SHA-1 digest.
# CTD data are subject to the CTD terms of use (https://ctdbase.org/about/legal.jsp).
#
# Run time: about 10 minutes, most of it in import_CTD(). Needs about 300 MB
# and the openssl package for the digest check.

library(ctdR)
library(data.table)

# A private cache, so the 2023 release does not replace the user's own.
work_dir <- file.path(tempdir(), "ctdR_yang2024")
dir.create(work_dir, showWarnings = FALSE)
options(ctdR.cache = file.path(work_dir, "cache"))

## 1. The CTD release of 23 May 2023, verified against the archive digest ----
ctd_file <- file.path(work_dir, "CTD_chem_gene_ixns_2023-05-23.csv.gz")
archive_capture_url <- paste0("https://web.archive.org/web/20230607020439id_/",
    "https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz")
if (!file.exists(ctd_file))
    utils::download.file(archive_capture_url, ctd_file, mode = "wb")
# The archive's index (CDX) records the digest of the capture it served.
archive_digest <- readLines(paste0("https://web.archive.org/cdx/search/cdx?url=",
    "ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz",
    "&from=20230607020439&to=20230607020439&fl=digest"))

#' Encode bytes in base32 (RFC 4648), the alphabet the archive prints.
#'
#' openssl returns the SHA-1 as raw bytes, the archive index as base32
#' text; this puts the local digest in the archive's form so the two
#' can be compared as strings.
#'
#' @param digest_bytes Raw vector, the bytes to encode (20 for a SHA-1,
#'   whose 160 bits split exactly into 32 base32 characters).
#' @return Character of length 1, the base32 encoding without padding.
base32 <- function(digest_bytes) {
    # most significant bit first within each byte
    digest_bits <- unlist(lapply(digest_bytes,
        function(one_byte) rev(as.integer(rawToBits(one_byte)))))
    # each group of five bits is one base32 character, value 0-31
    character_values <- vapply(
        split(digest_bits, ceiling(seq_along(digest_bits) / 5)),
        function(five_bits) sum(five_bits * 2^(4:0)), numeric(1))
    paste(c(LETTERS, 2:7)[character_values + 1], collapse = "")
}
local_digest <- base32(as.raw(openssl::sha1(file(ctd_file, raw = TRUE))))
stopifnot(identical(local_digest, archive_digest[1]))

## 2. Import every species: symbol sets pool species by gene symbol --------
import_CTD(ctd_file, organism = "all")
print(ctd_provenance())

## 3. Inputs published by the authors ---------------------------------------
yang_extdata_dir <- system.file("extdata", "yang2024", package = "ctdR")
dyslexia_genes <- fread(file.path(yang_extdata_dir,
    "yang2024_genes.tsv"))$GeneSymbol
published_results <- fread(file.path(yang_extdata_dir,
    "yang2024_replication.tsv"))

## 4. Same sets? Set size N and overlap k per chemical ----------------------
chemical_gene_sets <- as_genesets_CTD(id_type = "symbol")
ctdR_set_overlap <- data.table(ChemicalID = names(chemical_gene_sets),
    N_ctdR = lengths(lapply(chemical_gene_sets, unique)),
    k_ctdR = vapply(chemical_gene_sets,
        function(set_symbols) sum(unique(set_symbols) %in% toupper(dyslexia_genes)),
        0L))

## 5. The authors' test on ctdR's sets: hypergeometric, universe 19,000 ----
## Their p-values are reproduced with the 130 genes of their Table S2.

#' The authors' enrichment test: upper tail of the hypergeometric.
#'
#' Probability of drawing at least k_overlap genes of a chemical's set
#' when drawing the n_input dyslexia genes from a universe of
#' universe_size genes. Written out so the authors' settings stay
#' visible next to ctdR's own ORA.
#'
#' @param k_overlap Integer vector, dyslexia genes in the chemical's set.
#' @param N_set_size Integer vector, genes in the chemical's set.
#' @param n_input Integer, dyslexia genes in the query (130, Table S2).
#' @param universe_size Integer, genes in the universe (19,000, as the
#'   authors state).
#' @return Numeric vector of p-values; NA where k_overlap is 0, because
#'   the authors report only chemicals that share at least one gene.
hypergeometric_upper_tail <- function(k_overlap, N_set_size,
                                      n_input = 130, universe_size = 19000)
    ifelse(k_overlap > 0,
        phyper(k_overlap - 1, n_input, universe_size - n_input, N_set_size,
            lower.tail = FALSE),
        NA_real_)
replication_check <- merge(published_results[, .(Chemical, Class, ChemicalID,
    N_paper, k_paper, PValue_paper, PValueFDR_paper)], ctdR_set_overlap,
    by = "ChemicalID", all.x = TRUE)
replication_check[, PValue_sameTest_ctdR :=
    hypergeometric_upper_tail(k_ctdR, N_ctdR)]

## 6. Discovery: every chemical, ctdR's ORA --------------------------------
entrez_by_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    dyslexia_genes, "ENTREZID", "SYMBOL")
# Without alpha, ORA tests every row of the table, so all the dyslexia
# genes form the query; the pvalue column is a required placeholder, not a
# measurement.
ora_input <- data.frame(
    EntrezID = unname(entrez_by_symbol[!is.na(entrez_by_symbol)]),
    pvalue = 0.01)
protein_coding_symbols <- unique(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, keys = "protein-coding",
    keytype = "GENETYPE", columns = "SYMBOL")$SYMBOL)
all_human_symbols <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL")

#' Run ctdR's ORA on the dyslexia genes and rank every chemical.
#'
#' The three discovery analyses differ only in background and minimum set
#' size; one helper keeps everything else identical.
#'
#' @param universe Character vector, gene symbols of the background.
#' @param minimum_set_size Integer, smallest set tested, in genes.
#' @return data.table of the ORA result ordered by adjusted, then raw,
#'   p-value, with a Rank column (1 = most significant) and GeneSetSize,
#'   the number of the chemical's genes in the background.
run_discovery <- function(universe, minimum_set_size) {
    ora_result <- as.data.table(enrichment_CTD(ora_input, method = "ORA",
        universe = universe, minGSSize = minimum_set_size))
    setorder(ora_result, PValueAdjusted, PValue)
    ora_result[, `:=`(Rank = seq_len(.N),
        GeneSetSize = as.integer(sub("/.*", "", BackgroundRatio)))]
    ora_result
}
# main: human protein-coding background, the closest to the authors'
# 19,000; all_human: every human gene symbol; min10: sets of ten genes or
# more, the default minimum of pathway tools.
discovery_runs <- list(
    main = run_discovery(protein_coding_symbols, 2),
    all_human = run_discovery(all_human_symbols, 2),
    min10 = run_discovery(protein_coding_symbols, 10))
discovery <- discovery_runs$main
print(head(discovery[, .(ChemicalName, Count, FoldEnrichment, PValue,
    PValueAdjusted)], 30))

## 7. The tables of the vignette ---------------------------------------------
## Written to ./yang2024_outputs/; compare them with inst/extdata/yang2024/
## before replacing the bundled copies. Counts and p-values only, no CTD
## records. P-values keep three or four significant digits.
SIGNIFICANCE_THRESHOLD <- 0.05
output_dir <- file.path(getwd(), "yang2024_outputs")
dir.create(output_dir, showWarnings = FALSE)

#' Write one table of the vignette as TSV in the output directory.
#'
#' @param table_to_write data.frame or data.table.
#' @param file_name Character, file name with extension.
write_vignette_table <- function(table_to_write, file_name) {
    # no quotes: no field holds a tab, and an empty identifier stays empty
    fwrite(table_to_write, file.path(output_dir, file_name), sep = "\t",
        quote = FALSE)
}

## 7a. Same sets, same test: the paper's columns as published, in the
##     paper's order, next to the ones ctdR recomputes.
ctdR_row <- match(published_results$ChemicalID, replication_check$ChemicalID)
replication_table <- published_results[, .(Class, Chemical, ChemicalID,
    Mapping, N_paper, k_paper, PValue_paper, PValueFDR_paper)]
replication_table[, `:=`(
    N_ctdR = replication_check$N_ctdR[ctdR_row],
    k_ctdR = replication_check$k_ctdR[ctdR_row],
    PValue_sameTest_ctdR =
        signif(replication_check$PValue_sameTest_ctdR[ctdR_row], 4))]
write_vignette_table(replication_table, "yang2024_replication.tsv")

## 7b. Discovery: how many chemicals each analysis tests and finds.
write_vignette_table(rbindlist(lapply(names(discovery_runs),
    function(analysis_name) {
        ranking <- discovery_runs[[analysis_name]]
        data.table(Analysis = analysis_name, Tested = nrow(ranking),
            Significant = sum(ranking$PValueAdjusted < SIGNIFICANCE_THRESHOLD),
            NominalP05 = sum(ranking$PValue < SIGNIFICANCE_THRESHOLD))
    })), "yang2024_discovery_summary.tsv")

## 7c. The top 30 of the main analysis, flagged when the chemical is in the
##     authors' list.
paper_chemical_ids <- unique(published_results[!is.na(ChemicalID) &
    nzchar(ChemicalID), ChemicalID])
write_vignette_table(discovery[1:30, .(Rank, ChemicalID, ChemicalName, Count,
    GeneSetSize, FoldEnrichment = round(FoldEnrichment, 2),
    PValue = signif(PValue, 3), PValueAdjusted = signif(PValueAdjusted, 3),
    InYangList = ChemicalID %in% paper_chemical_ids)],
    "yang2024_discovery_top30.tsv")

## 7d. Where the chemicals the authors found significant fall, in each
##     analysis. A chemical that appears in two classes counts once.
yang_significant <- unique(published_results[!is.na(PValueFDR_paper) &
    PValueFDR_paper < SIGNIFICANCE_THRESHOLD, .(ChemicalID, Chemical, Class)],
    by = "ChemicalID")
yang_positions <- rbindlist(lapply(names(discovery_runs),
    function(analysis_name) {
        ranking <- discovery_runs[[analysis_name]]
        positions <- ranking[match(yang_significant$ChemicalID,
            ranking$ChemicalID)]
        data.table(yang_significant, Analysis = analysis_name,
            Rank = positions$Rank, Count = positions$Count,
            GeneSetSize = positions$GeneSetSize,
            PValue = signif(positions$PValue, 3),
            PValueAdjusted = signif(positions$PValueAdjusted, 3),
            # blank where the chemical was not tested
            Tested = ifelse(is.na(positions$Rank), NA_integer_,
                nrow(ranking)))
    }))
setorder(yang_positions, Analysis, Rank, na.last = TRUE)
write_vignette_table(yang_positions, "yang2024_discovery_yang_chemicals.tsv")

## 7e. Yang et al. against ctdR, chemical by chemical. One row per entry of
##     the paper's tables: the authors' call next to ctdR's call in the main
##     discovery analysis, both at 0.05, twice. Unadjusted: the paper's
##     p-value against ctdR's raw p-value. Adjusted: the paper's FDR (within
##     each class, as published) against ctdR's Benjamini-Hochberg p-value
##     (over every chemical tested).

#' Compare the authors' call with ctdR's for each chemical.
#'
#' The same rule is applied with and without multiple-testing correction,
#' so the two cross-counts differ only in the p-values they read.
#'
#' @param paper_pvalue Numeric vector, the authors' p-value or FDR per
#'   chemical; NA where the paper gives none, which happens only when no
#'   dyslexia gene is in the set (k = 0), so NA counts as not significant.
#' @param ctdR_pvalue Numeric vector, ctdR's p-value for the same chemicals;
#'   NA where ctdR did not test the chemical.
#' @return Character vector, one of "significant in both",
#'   "significant only in Yang et al.", "significant only in ctdR",
#'   "significant in neither" or "not tested by ctdR". A chemical ctdR did
#'   not test is reported as such, not as non-significant.
compare_calls <- function(paper_pvalue, ctdR_pvalue) {
    significant_in_paper <- !is.na(paper_pvalue) &
        paper_pvalue < SIGNIFICANCE_THRESHOLD
    significant_in_ctdR <- ctdR_pvalue < SIGNIFICANCE_THRESHOLD
    ifelse(is.na(ctdR_pvalue), "not tested by ctdR",
    ifelse(significant_in_paper & significant_in_ctdR, "significant in both",
    ifelse(significant_in_paper, "significant only in Yang et al.",
    ifelse(significant_in_ctdR, "significant only in ctdR",
        "significant in neither"))))
}

discovery_row <- match(replication_check$ChemicalID, discovery$ChemicalID)
yang_vs_ctdR <- replication_check[, .(Class, Chemical, ChemicalID,
    k_paper, N_paper, PValue_paper, PValueFDR_paper, PValue_sameTest_ctdR)]
yang_vs_ctdR[, `:=`(
    Rank_discovery = discovery$Rank[discovery_row],
    Tested_discovery = nrow(discovery),
    PValue_discovery = discovery$PValue[discovery_row],
    PValueAdjusted_discovery = discovery$PValueAdjusted[discovery_row])]
yang_vs_ctdR[, `:=`(
    Outcome_unadjusted = compare_calls(PValue_paper, PValue_discovery),
    Outcome_adjusted = compare_calls(PValueFDR_paper,
        PValueAdjusted_discovery))]
setorder(yang_vs_ctdR, PValueFDR_paper, Rank_discovery, na.last = TRUE)
yang_vs_ctdR[, `:=`(
    PValue_sameTest_ctdR = signif(PValue_sameTest_ctdR, 4),
    PValue_discovery = signif(PValue_discovery, 4),
    PValueAdjusted_discovery = signif(PValueAdjusted_discovery, 4))]
write_vignette_table(yang_vs_ctdR, "yang2024_comparison.tsv")
print(table(Unadjusted = yang_vs_ctdR$Outcome_unadjusted,
    Adjusted = yang_vs_ctdR$Outcome_adjusted))
message("Tables written to ", output_dir)
