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
discovery <- enrichment_CTD(ora_input, method = "ORA",
    universe = protein_coding_symbols, minGSSize = 2)
head(discovery[, c("ChemicalName", "Count", "FoldEnrichment",
    "PValue", "PValueAdjusted")], 30)

## 7. Yang et al. against ctdR, chemical by chemical --------------------------
## One row per chemical of the paper's tables: the authors' call next to
## ctdR's call in discovery, both at 0.05, twice. Unadjusted: the paper's
## p-value against ctdR's raw p-value. Adjusted: the paper's FDR (within each
## class, as published) against ctdR's Benjamini-Hochberg p-value (over every
## chemical tested). The vignette shows the two cross-counts and the rows
## where the calls disagree.
SIGNIFICANCE_THRESHOLD <- 0.05
output_dir <- file.path(getwd(), "yang2024_outputs")
dir.create(output_dir, showWarnings = FALSE)

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

discovery_by_pvalue <- discovery[order(discovery$PValue), ]
discovery_row <- match(replication_check$ChemicalID,
    discovery_by_pvalue$ChemicalID)
yang_vs_ctdR <- replication_check[, .(Class, Chemical, ChemicalID,
    k_paper, N_paper, PValue_paper, PValueFDR_paper, PValue_sameTest_ctdR)]
yang_vs_ctdR[, `:=`(
    Rank_discovery = discovery_row,
    Tested_discovery = nrow(discovery_by_pvalue),
    PValue_discovery = discovery_by_pvalue$PValue[discovery_row],
    PValueAdjusted_discovery =
        discovery_by_pvalue$PValueAdjusted[discovery_row])]
yang_vs_ctdR[, `:=`(
    Outcome_unadjusted = compare_calls(PValue_paper, PValue_discovery),
    Outcome_adjusted = compare_calls(PValueFDR_paper,
        PValueAdjusted_discovery))]
setorder(yang_vs_ctdR, PValueFDR_paper, Rank_discovery, na.last = TRUE)
# four significant digits, as in the other bundled tables
yang_vs_ctdR[, `:=`(
    PValue_sameTest_ctdR = signif(PValue_sameTest_ctdR, 4),
    PValue_discovery = signif(PValue_discovery, 4),
    PValueAdjusted_discovery = signif(PValueAdjusted_discovery, 4))]
fwrite(yang_vs_ctdR, file.path(output_dir, "yang2024_comparison.tsv"),
    sep = "\t")
print(table(Unadjusted = yang_vs_ctdR$Outcome_unadjusted,
    Adjusted = yang_vs_ctdR$Outcome_adjusted))
message("Comparison written to ", file.path(output_dir,
    "yang2024_comparison.tsv"))
