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
## Required packages: org.Hs.eg.db, AnnotationDbi, SummarizedExperiment
##          (already in Imports of ctdR)

suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(SummarizedExperiment)
})

GEO_URL <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311566/",
    "suppl/GSE311566_PBMCs_Female_normalized_counts.txt.gz"
)
## Number of genes kept by variance across the 7 samples: enough for the
## vignette to show every method, small enough to keep the package light.
N_TOP_VARIANCE_GENES <- 1500L

## 1. Download to tempdir (kept in cache for reruns) ------------------------
counts_gz_path <- file.path(tempdir(),
    "GSE311566_Female_normalized_counts.txt.gz")
if (!file.exists(counts_gz_path)) {
    message("Downloading ", GEO_URL)
    utils::download.file(GEO_URL, counts_gz_path, mode = "wb", quiet = FALSE)
}

## 2. Read counts (ENSG rownames, NAME + DESCRIPTION + 20 sample columns) ---
normalised_counts <- utils::read.table(
    gzfile(counts_gz_path),
    header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
)
normalised_counts$DESCRIPTION <- NULL  # all NA in this dataset
message("Loaded ", nrow(normalised_counts), " genes x ",
        ncol(normalised_counts), " samples")

## 3. Keep only Ctrl_F_* and DEX_F_* columns -------------------------------
## The file also holds the PFAS arms; the example uses vehicle and Dex only.
dex_vs_dmso_columns <- grep("^(Ctrl|DEX)_F_", colnames(normalised_counts),
    value = TRUE)
stopifnot(length(dex_vs_dmso_columns) >= 6)  # need at least 3+3 for t-test
normalised_counts <- normalised_counts[, dex_vs_dmso_columns]
treatment_group <- factor(
    ifelse(grepl("^DEX_", dex_vs_dmso_columns), "Dex", "DMSO"),
    levels = c("DMSO", "Dex")
)
samples_per_group <- table(treatment_group)
message("Kept ", length(dex_vs_dmso_columns), " samples: ",
        paste(samples_per_group, names(samples_per_group),
              sep = " ", collapse = ", "))

## 4. Map ENSG -> Entrez ----------------------------------------------------
## Strip the version suffix (ENSG00000141510.17 -> ENSG00000141510),
## which org.Hs.eg.db does not carry.
ensembl_ids <- sub("\\..*$", "", rownames(normalised_counts))
entrez_by_ensembl <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids, keytype = "ENSEMBL", column = "ENTREZID",
    multiVals = "first"
)
## CTD indexes genes by Entrez ID: drop genes without one, and keep the
## first of several Ensembl IDs mapping to the same Entrez ID.
has_unique_entrez <- !is.na(entrez_by_ensembl) &
    !duplicated(entrez_by_ensembl)
normalised_counts <- normalised_counts[has_unique_entrez, ]
rownames(normalised_counts) <- entrez_by_ensembl[has_unique_entrez]
message("After Entrez mapping: ", nrow(normalised_counts), " genes")

## 5. log2(x+1) transform on already-normalised counts ---------------------
log2_expression <- log2(as.matrix(normalised_counts) + 1)

## 6. Force-include CTD-sample target Entrez IDs (so the example actually
##    produces non-empty ORA/CAMERA results on the bundled chemical set) ---
CTD_SAMPLE_ENTREZ <- as.character(c(
    1956, 207, 2099, 3553, 3569, 3845, 4609, 4790,
    5290, 595, 596, 672, 6774, 7124, 7157, 7422, 836
))
ctd_sample_genes_present <- intersect(CTD_SAMPLE_ENTREZ,
    rownames(log2_expression))
message("Force-keeping ", length(ctd_sample_genes_present), " of ",
        length(CTD_SAMPLE_ENTREZ), " CTD-sample target genes")

## 7. Top-variance filter ---------------------------------------------------
gene_variance <- apply(log2_expression, 1, stats::var)
top_variance_genes <- names(sort(gene_variance, decreasing = TRUE))[
    seq_len(min(N_TOP_VARIANCE_GENES, length(gene_variance)))]
subset_genes <- union(ctd_sample_genes_present, top_variance_genes)
log2_expression <- log2_expression[subset_genes, , drop = FALSE]
message("Final subset: ", nrow(log2_expression), " genes x ",
        ncol(log2_expression), " samples")

## 8. Save -------------------------------------------------------------------
## Stored as a SummarizedExperiment: the assay and the sample annotation stay
## in one object, so subsetting or reordering samples cannot silently
## desynchronise them from the group labels used to build the design matrix.
## The column is named "group": the vignette and the man pages read gse$group.
sample_annotation <- DataFrame(
    sample = colnames(log2_expression),
    group = treatment_group,
    row.names = colnames(log2_expression)
)
gse_subset <- SummarizedExperiment(
    assays = list(logcounts = log2_expression),
    colData = sample_annotation,
    metadata = list(
        source = "GEO GSE311566",
        url = GEO_URL,
        subset = paste0(
            "Female donors, Dexamethasone vs vehicle (Ctrl/DMSO); ",
            nrow(log2_expression), " genes x ", ncol(log2_expression),
            " samples"
        ),
        assay_units = "log2(normalised count + 1)",
        generated_by = "inst/scripts/make_gse311566_subset.R"
    )
)
output_path <- file.path("inst", "extdata", "GSE311566_subset.rds")
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(gse_subset, file = output_path, compress = "xz")

message("Wrote ", output_path, " (",
        round(file.info(output_path)$size / 1024, 1), " KB)")
