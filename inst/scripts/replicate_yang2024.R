# Replicate Yang et al. (2024) with ctdR, then screen the whole CTD catalogue
# with the same genes. Produces the tables in inst/extdata/yang2024/ that the
# vignette "Replicating a published CTD analysis" displays.
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

work <- file.path(tempdir(), "ctdR_yang2024")
dir.create(work, showWarnings = FALSE)
options(ctdR.cache = file.path(work, "cache"))

## 1. The CTD release of 23 May 2023, verified against the archive digest ----
ctd_file <- file.path(work, "CTD_chem_gene_ixns_2023-05-23.csv.gz")
capture <- paste0("https://web.archive.org/web/20230607020439id_/",
    "https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz")
if (!file.exists(ctd_file))
    utils::download.file(capture, ctd_file, mode = "wb")
cdx <- readLines(paste0("https://web.archive.org/cdx/search/cdx?url=",
    "ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz",
    "&from=20230607020439&to=20230607020439&fl=digest"))
# The archive prints SHA-1 digests in base32 (RFC 4648).
base32 <- function(raw) {
    # most significant bit first within each byte
    bits <- unlist(lapply(raw, function(b) rev(as.integer(rawToBits(b)))))
    idx <- vapply(split(bits, ceiling(seq_along(bits) / 5)),
        function(b) sum(b * 2^(4:0)), numeric(1))
    paste(c(LETTERS, 2:7)[idx + 1], collapse = "")
}
sha1 <- base32(as.raw(openssl::sha1(file(ctd_file, raw = TRUE))))
stopifnot(identical(sha1, cdx[1]))

## 2. Import every species: symbol sets pool species by gene symbol --------
import_CTD(ctd_file, organism = "all")
print(ctd_provenance())

## 3. Inputs published by the authors ---------------------------------------
ext <- system.file("extdata", "yang2024", package = "ctdR")
genes <- fread(file.path(ext, "yang2024_genes.tsv"))$GeneSymbol
rep <- fread(file.path(ext, "yang2024_replication.tsv"))

## 4. Same sets? Set size N and overlap k per chemical ----------------------
sets <- as_genesets_CTD(id_type = "symbol")
size <- data.table(ChemicalID = names(sets),
    N_ctdR = lengths(lapply(sets, unique)),
    k_ctdR = vapply(sets, function(g) sum(unique(g) %in% toupper(genes)), 0L))

## 5. The authors' test on ctdR's sets: hypergeometric, universe 19,000 ----
## Their p-values are reproduced with the 130 genes of their Table S2.
hyper <- function(k, N, n = 130, M = 19000)
    ifelse(k > 0, phyper(k - 1, n, M - n, N, lower.tail = FALSE), NA_real_)
check <- merge(rep[, .(Chemical, Class, ChemicalID, N_paper, k_paper,
    PValue_paper, PValueFDR_paper)], size, by = "ChemicalID", all.x = TRUE)
check[, PValue_sameTest_ctdR := hyper(k_ctdR, N_ctdR)]

## 6. Discovery: every chemical, ctdR's ORA --------------------------------
entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, genes,
    "ENTREZID", "SYMBOL")
x <- data.frame(EntrezID = unname(entrez[!is.na(entrez)]), pvalue = 0.01)
coding <- unique(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = "protein-coding", keytype = "GENETYPE", columns = "SYMBOL")$SYMBOL)
discovery <- enrichment_CTD(x, method = "ORA", universe = coding, minGSSize = 2)
head(discovery[, c("ChemicalName", "Count", "FoldEnrichment",
    "PValue", "PValueAdjusted")], 30)
