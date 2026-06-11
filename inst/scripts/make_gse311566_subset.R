## Generate inst/extdata/GSE311566_subset.rds
##
## Source : Comparing the effects of PFAS and dexamethasone on human PBMCs
##          GEO accession GSE311566
##          https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311566
##
## Subset : Female donors, Dex vs vehicle (Ctrl/DMSO) contrast only
##          7 samples (4 Ctrl + 3 Dex), top-variance genes mapped to Entrez
##
## Purpose: bundled real-data example for ctdR vignette and man pages.
##          NOT for scientific re-analysis -- subset is filtered for size.
##
## Licence: NCBI GEO public submission. Users of this subset should cite the
##          original GSE311566 series. See inst/extdata/README.md for details.
##
## To regenerate the bundled .rds, run from package root:
##   Rscript inst/scripts/make_gse311566_subset.R
##
## Required packages: org.Hs.eg.db, AnnotationDbi (already in Imports of ctdR)

suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})

GEO_URL <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311566/",
    "suppl/GSE311566_PBMCs_Female_normalized_counts.txt.gz"
)
N_TOP_VAR <- 1500L

## 1. Download to tempdir (kept in cache for reruns) ------------------------
tmp_gz <- file.path(tempdir(), "GSE311566_Female_normalized_counts.txt.gz")
if (!file.exists(tmp_gz)) {
    message("Downloading ", GEO_URL)
    utils::download.file(GEO_URL, tmp_gz, mode = "wb", quiet = FALSE)
}

## 2. Read counts (ENSG rownames, NAME + DESCRIPTION + 20 sample columns) ---
counts <- utils::read.table(
    gzfile(tmp_gz),
    header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
)
counts$DESCRIPTION <- NULL  # all NA in this dataset
message("Loaded ", nrow(counts), " genes x ", ncol(counts), " samples")

## 3. Keep only Ctrl_F_* and DEX_F_* columns -------------------------------
keep <- grep("^(Ctrl|DEX)_F_", colnames(counts), value = TRUE)
stopifnot(length(keep) >= 6)  # need at least 3+3 for t-test
counts <- counts[, keep]
group <- factor(
    ifelse(grepl("^DEX_", keep), "Dex", "DMSO"),
    levels = c("DMSO", "Dex")
)
message("Kept ", length(keep), " samples: ",
        paste(table(group), names(table(group)), sep = " ", collapse = ", "))

## 4. Map ENSG -> Entrez ----------------------------------------------------
ensg <- sub("\\..*$", "", rownames(counts))  # strip version if any
entrez_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensg, keytype = "ENSEMBL", column = "ENTREZID",
    multiVals = "first"
)
keep_g <- !is.na(entrez_map) & !duplicated(entrez_map)
counts <- counts[keep_g, ]
rownames(counts) <- entrez_map[keep_g]
message("After Entrez mapping: ", nrow(counts), " genes")

## 5. log2(x+1) transform on already-normalised counts ---------------------
expr <- log2(as.matrix(counts) + 1)

## 6. Force-include CTD-sample target Entrez IDs (so the example actually
##    produces non-empty ORA/CAMERA results on the bundled chemical set) ---
CTD_SAMPLE_ENTREZ <- as.character(c(
    1956, 207, 2099, 3553, 3569, 3845, 4609, 4790,
    5290, 595, 596, 672, 6774, 7124, 7157, 7422, 836
))
must_keep <- intersect(CTD_SAMPLE_ENTREZ, rownames(expr))
message("Force-keeping ", length(must_keep), " of ",
        length(CTD_SAMPLE_ENTREZ), " CTD-sample target genes")

## 7. Top-variance filter ---------------------------------------------------
vars <- apply(expr, 1, stats::var)
top_var <- names(sort(vars, decreasing = TRUE))[seq_len(min(N_TOP_VAR, length(vars)))]
final_genes <- union(must_keep, top_var)
expr <- expr[final_genes, , drop = FALSE]
message("Final subset: ", nrow(expr), " genes x ", ncol(expr), " samples")

## 8. Save -----------------------------------------------------------------
coldata <- data.frame(
    sample = colnames(expr),
    group = group,
    row.names = colnames(expr),
    stringsAsFactors = FALSE
)
out <- list(expr = expr, coldata = coldata)
out_path <- file.path("inst", "extdata", "GSE311566_subset.rds")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file = out_path, compress = "xz")

message("Wrote ", out_path, " (",
        round(file.info(out_path)$size / 1024, 1), " KB)")
