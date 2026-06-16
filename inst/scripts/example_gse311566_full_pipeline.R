## ============================================================
## ctdR end-to-end example -- full GSE311566 Dex-vs-DMSO pipeline
## ============================================================
##
## Purpose
##   A runnable, self-contained reference pipeline that exercises the
##   *full* ctdR workflow on a real RNA-seq dataset, with explicit
##   a-priori statistical power analysis and declared significance
##   thresholds. Intended as the canonical example linked from the
##   paper and the GitHub README.
##
## Differences from the vignette
##   The package vignette uses a *subset* (1,500 top-variance genes,
##   7 samples) of GSE311566 packaged inside `inst/extdata/`. The
##   vignette reports raw p-values *without* applying alpha cutoffs
##   because the subset is too small to support a significance claim.
##   This script instead uses the *full* normalised-count matrix
##   from GEO and applies declared alpha thresholds throughout, as
##   one would in a production analysis.
##
## Scope and limits (Steps A->E)
##   This pipeline is an *illustration* of the ctdR API on real data.
##   To keep it self-contained and readable, Steps A->E apply only
##   the bare minimum statistical machinery -- limma `lmFit + eBayes`
##   on log2(normalised counts + 1), followed by BH multiple-testing
##   correction. None of the corrections a real biomarker-discovery
##   study would apply are performed here, in particular:
##     * no batch-effect correction (e.g. ComBat / SVA / RUVseq);
##     * no adjustment for covariates (donor, age, library prep);
##     * no per-sample quality control or outlier exclusion;
##     * no expression-level filtering of low-count genes;
##     * no robustness checks (resampling, permutation null).
##   A real analysis would account for these. Do not read the numbers
##   produced here as a scientific claim about Dex response in PBMCs
##   -- read them as "what ctdR returns when handed a plausibly
##   prepared DE result on a public dataset".
##
## Dataset
##   GEO accession GSE311566 -- "Comparing the effects of PFAS and
##   dexamethasone on human PBMCs" (Female donors). We restrict to
##   the Dex vs vehicle (Ctrl/DMSO) contrast.
##   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311566
##
## Prerequisites
##   1. Internet access to download the GEO supplementary file
##      (~5 MB compressed).
##   2. A populated ctdR cache. You must have run import_CTD() once
##      against the real CTD chemical-gene interactions file:
##        download CTD_chem_gene_ixns.csv.gz from
##        https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz
##        gunzip and: ctdR::import_CTD("path/to/CTD_chem_gene_ixns.csv")
##      The bundled inst/extdata/CTD_chem_gene_ixns_sample.csv is
##      *not* used by this script (it is a 10-chemical toy file that
##      cannot support a significance claim).
##
## Run
##   Rscript inst/scripts/example_gse311566_full_pipeline.R
##
## Outputs
##   Written to ./example_outputs/ (created if missing):
##     - de_full.tsv                  per-gene DE table
##     - {ora,gsea,camera}_full.tsv   all chemicals returned (for inspection)
##     - {ora,gsea,camera}_significant.tsv  subset passing alpha_fdr_chemical
##     - gsva_scores.tsv              GSVA chemical x sample matrix
##     - plot_{ora,gsea,camera,gsva}.png  one bar/heatmap per method
##     - power_report.txt             a-priori power analysis summary
##     - top_chemicals_summary.tsv    top-N hits per method by raw p-value
##     - expected_hit_ranks.tsv       where Dexamethasone (D003907) lands
##     - sessionInfo.txt              R session info for reproducibility
##
## Licensing
##   This script downloads public GEO data and reads cached CTD data
##   the user is responsible for obtaining. ctdR does not bundle or
##   redistribute CTD data. See https://ctdbase.org/about/legal.jsp.
## ============================================================

suppressPackageStartupMessages({
    library(ctdR)
    library(limma)
    library(ggplot2)
    library(rappdirs)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

## ---- Config ------------------------------------------------
## All knobs live here. Edit and rerun to explore.

CONFIG <- list(
    geo_url = paste0(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311566/",
        "suppl/GSE311566_PBMCs_Female_normalized_counts.txt.gz"
    ),

    ## Significance thresholds (declared a priori).
    alpha_nominal_de   = 0.05,  # per-gene t-test nominal alpha
    alpha_fdr_de       = 0.05,  # BH-adjusted alpha for DE calls
    alpha_fdr_chemical = 0.05,  # BH-adjusted alpha for chemical hits

    ## Effect-size target for the a-priori power calculation.
    ## log2FC = 1 (a doubling) with residual sd ~ 0.5 on log2 scale
    ## is conservative for known Dex-responsive genes in PBMCs.
    power_log2fc_target = 1.0,
    power_sd_target     = 0.5,

    ## Multiple-testing correction passed through to enrichment_CTD().
    p_adjust_method = "BH",

    out_dir = file.path(getwd(), "example_outputs")
)

dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)

log_step <- function(...) {
    message(sprintf("[%s] %s",
        format(Sys.time(), "%H:%M:%S"),
        paste0(..., collapse = "")
    ))
}

## ---- Step A -- a-priori power analysis ---------------------
##
## With n1 = 4 DMSO, n2 = 3 Dex, what is the per-gene power of a
## two-sample t-test to detect log2FC = 1 at sd = 0.5 (alpha = 0.05,
## two-sided)? power.t.test() assumes equal-n; we use the harmonic
## mean as a slightly conservative effective n per group.

n1 <- 4L; n2 <- 3L
n_harm <- 2 / (1 / n1 + 1 / n2)

pwr <- stats::power.t.test(
    n     = n_harm,
    delta = CONFIG$power_log2fc_target,
    sd    = CONFIG$power_sd_target,
    sig.level = CONFIG$alpha_nominal_de,
    type  = "two.sample",
    alternative = "two.sided"
)

power_report <- c(
    sprintf("ctdR full-pipeline example -- a-priori power analysis"),
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf(""),
    sprintf("Design : n1 = %d DMSO, n2 = %d Dex (harmonic n = %.2f)", n1, n2, n_harm),
    sprintf("Target : log2FC = %.2f, residual sd = %.2f",
        CONFIG$power_log2fc_target, CONFIG$power_sd_target),
    sprintf("Alpha  : %.3f (two-sided, per-gene)", CONFIG$alpha_nominal_de),
    sprintf("Power  : %.3f", pwr$power),
    sprintf(""),
    sprintf("Interpretation: ~%.0f%% probability of declaring a true ", pwr$power * 100),
    sprintf("log2FC=%.1f gene significant before multiple-testing ",
        CONFIG$power_log2fc_target),
    sprintf("correction. After BH at FDR <= %.2f the realized ",
        CONFIG$alpha_fdr_de),
    sprintf("power per gene is lower; enrichment methods recover ranking "),
    sprintf("signal even when most individual genes do not pass FDR.")
)
writeLines(power_report, file.path(CONFIG$out_dir, "power_report.txt"))
log_step("Step A -- power analysis: power = ", sprintf("%.3f", pwr$power))

## ---- Step B -- download + load full GSE311566 --------------

tmp_gz <- file.path(tempdir(),
    "GSE311566_PBMCs_Female_normalized_counts.txt.gz")
if (!file.exists(tmp_gz)) {
    log_step("Step B -- downloading ", CONFIG$geo_url)
    utils::download.file(CONFIG$geo_url, tmp_gz, mode = "wb", quiet = TRUE)
}

counts <- utils::read.table(
    gzfile(tmp_gz),
    header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
)
counts$DESCRIPTION <- NULL

keep_cols <- grep("^(Ctrl|DEX)_F_", colnames(counts), value = TRUE)
stopifnot(length(keep_cols) >= 6)
counts <- counts[, keep_cols]

group <- factor(
    ifelse(grepl("^DEX_", keep_cols), "Dex", "DMSO"),
    levels = c("DMSO", "Dex")
)
log_step("Step B -- loaded ", nrow(counts), " genes x ",
    ncol(counts), " samples (",
    paste(table(group), names(table(group)), sep = " ", collapse = ", "),
    ")")

## ---- Step C -- map ENSG -> Entrez --------------------------

ensg <- sub("\\..*$", "", rownames(counts))
entrez_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensg, keytype = "ENSEMBL", column = "ENTREZID",
    multiVals = "first"
)
keep_g <- !is.na(entrez_map) & !duplicated(entrez_map)
counts <- counts[keep_g, ]
rownames(counts) <- entrez_map[keep_g]
log_step("Step C -- mapped to Entrez: ", nrow(counts), " genes")

expr <- log2(as.matrix(counts) + 1)

## ---- Step D -- differential expression with limma ----------
##
## Data are already library-normalised counts. We fit lmFit + eBayes
## on the log2(x+1) matrix. For raw-count data prefer limma-voom or
## DESeq2; here we follow the GEO supplement's already-normalised
## form to keep the example single-file.

design <- model.matrix(~ group)
fit <- limma::lmFit(expr, design)
fit <- limma::eBayes(fit)
tt <- limma::topTable(fit, coef = "groupDex",
    number = Inf, sort.by = "P")

de <- data.frame(
    EntrezID = rownames(tt),
    log2FC   = tt$logFC,
    pvalue   = tt$P.Value,
    padj     = stats::p.adjust(tt$P.Value, method = CONFIG$p_adjust_method),
    row.names = NULL,
    stringsAsFactors = FALSE
)
utils::write.table(de, file.path(CONFIG$out_dir, "de_full.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

n_nominal <- sum(de$pvalue < CONFIG$alpha_nominal_de)
n_fdr     <- sum(de$padj   < CONFIG$alpha_fdr_de)
log_step(sprintf(
    "Step D -- DE: %d genes pvalue < %.2f, %d genes BH-padj < %.2f",
    n_nominal, CONFIG$alpha_nominal_de, n_fdr, CONFIG$alpha_fdr_de
))

## ---- Step E -- verify CTD cache is populated with real data ---

cache_dir <- rappdirs::user_cache_dir("ctdR")
chem_rda <- file.path(cache_dir, "chemicals.rda")
if (!file.exists(file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda")) ||
    !file.exists(chem_rda)) {
    stop(
        "CTD data not found in ", cache_dir, ".\n",
        "Download CTD_chem_gene_ixns.csv.gz from\n",
        "  https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz\n",
        "gunzip and run:\n",
        "  ctdR::import_CTD(\"path/to/CTD_chem_gene_ixns.csv\")\n",
        "This script intentionally refuses the bundled toy sample.",
        call. = FALSE
    )
}

## The bundled toy sample has 10 chemicals; the real CTD release has
## tens of thousands. Refuse anything that looks like the toy cache,
## so we never silently report "0 chemicals at FDR < 0.05" coming
## from a 10-chemical universe.
MIN_REAL_CHEMICALS <- 1000L
e <- new.env(parent = emptyenv())
load(chem_rda, envir = e)
n_cached <- nrow(e$chemicals)
if (n_cached < MIN_REAL_CHEMICALS) {
    stop(
        "CTD cache at ", cache_dir, " holds only ", n_cached,
        " chemicals -- this looks like the bundled toy sample.\n",
        "This script requires the real CTD release. Download\n",
        "  https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz\n",
        "gunzip and re-import:\n",
        "  ctdR::import_CTD(\"path/to/CTD_chem_gene_ixns.csv\")\n",
        "before re-running.",
        call. = FALSE
    )
}
log_step("Step E -- CTD cache OK at ", cache_dir,
    " (", n_cached, " chemicals)")

## ---- Step F -- enrichment with all four methods ------------

apply_alpha <- function(df, alpha) {
    df[!is.na(df$PValueAdjusted) & df$PValueAdjusted < alpha, , drop = FALSE]
}

## Flatten list-columns (e.g. fgsea's LeadingEdge) to comma-separated
## strings so the result can be written to a TSV.
flatten_list_cols <- function(df) {
    for (col in colnames(df)) {
        if (is.list(df[[col]])) {
            df[[col]] <- vapply(df[[col]],
                function(x) paste(x, collapse = ","),
                character(1))
        }
    }
    df
}

write_method <- function(df, name) {
    utils::write.table(flatten_list_cols(df),
        file.path(CONFIG$out_dir, paste0(name, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE)
}

## ORA: significant genes only (BH-FDR < alpha_fdr_de)
sig_de <- de[de$padj < CONFIG$alpha_fdr_de,
    c("EntrezID", "pvalue"), drop = FALSE]
if (nrow(sig_de) == 0) {
    log_step("Step F -- no genes pass FDR for ORA; ",
        "falling back to nominal p < ",
        CONFIG$alpha_nominal_de, " as input set")
    sig_de <- de[de$pvalue < CONFIG$alpha_nominal_de,
        c("EntrezID", "pvalue"), drop = FALSE]
}
ora_full <- enrichment_CTD(sig_de, method = "ORA",
    pAdjustMethod = CONFIG$p_adjust_method)
ora_sig <- apply_alpha(ora_full, CONFIG$alpha_fdr_chemical)
write_method(ora_full, "ora_full")
write_method(ora_sig, "ora_significant")
log_step("Step F -- ORA: ", nrow(ora_sig),
    " chemicals at FDR < ", CONFIG$alpha_fdr_chemical,
    " (", nrow(ora_full), " total returned)")

## GSEA: full ranked list
gsea_full <- enrichment_CTD(de[, c("EntrezID", "pvalue")],
    method = "GSEA",
    pAdjustMethod = CONFIG$p_adjust_method)
gsea_sig <- apply_alpha(gsea_full, CONFIG$alpha_fdr_chemical)
write_method(gsea_full, "gsea_full")
write_method(gsea_sig, "gsea_significant")
log_step("Step F -- GSEA: ", nrow(gsea_sig),
    " chemicals at FDR < ", CONFIG$alpha_fdr_chemical,
    " (", nrow(gsea_full), " total returned)")

## CAMERA: expression matrix + design + contrast (col 2 = Dex vs DMSO)
camera_full <- enrichment_CTD(expr, method = "CAMERA",
    design = design, contrast = 2,
    pAdjustMethod = CONFIG$p_adjust_method)
camera_sig <- apply_alpha(camera_full, CONFIG$alpha_fdr_chemical)
write_method(camera_full, "camera_full")
write_method(camera_sig, "camera_significant")
log_step("Step F -- CAMERA: ", nrow(camera_sig),
    " chemicals at FDR < ", CONFIG$alpha_fdr_chemical,
    " (", nrow(camera_full), " total returned)")

## GSVA: per-sample scores (no alpha; descriptive output)
gsva_scores <- enrichment_CTD(expr, method = "GSVA")
utils::write.table(
    data.frame(ChemicalID = rownames(gsva_scores), gsva_scores,
        check.names = FALSE),
    file.path(CONFIG$out_dir, "gsva_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
log_step("Step F -- GSVA: ", nrow(gsva_scores), " chemicals x ",
    ncol(gsva_scores), " samples")

## ---- Step F.5 -- top-N summary (printed + saved) -----------

TOP_N <- 5L
top_summary <- function(df, method) {
    if (!nrow(df)) return(NULL)
    df <- df[order(df$PValue), , drop = FALSE]
    data.frame(
        Method = method,
        Rank = seq_len(min(TOP_N, nrow(df))),
        ChemicalID = df$ChemicalID[seq_len(min(TOP_N, nrow(df)))],
        ChemicalName = df$ChemicalName[seq_len(min(TOP_N, nrow(df)))],
        PValue = df$PValue[seq_len(min(TOP_N, nrow(df)))],
        PValueAdjusted = df$PValueAdjusted[seq_len(min(TOP_N, nrow(df)))],
        stringsAsFactors = FALSE
    )
}
top_all <- do.call(rbind, list(
    top_summary(ora_full, "ORA"),
    top_summary(gsea_full, "GSEA"),
    top_summary(camera_full, "CAMERA")
))
utils::write.table(top_all,
    file.path(CONFIG$out_dir, "top_chemicals_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
message("\n--- Top ", TOP_N, " chemicals per method (by raw p-value) ---")
print(top_all, row.names = FALSE)
message("")

## ---- Step F.6 -- expected-hit ranking ----------------------
##
## The treatment in this contrast is Dexamethasone (CTD: D003907).
##
## On the bundled vignette subset the gene universe is built to
## *include* Dex's CTD targets, so Dex appears top-3 by GSEA. On the
## *unfiltered* full dataset against the *full* CTD release (~11k
## chemicals), Dex must compete with every co-regulated compound.
##
## Method differences observed here are characteristics of the
## methods, not limitations of ctdR. The four engines exposed by
## enrichment_CTD() trade sensitivity for specificity differently:
##   - GSEA uses the *full ranking* and recovers Dex even when most
##     individual DE genes do not pass FDR -- expected to be the
##     most sensitive on small-n studies.
##   - CAMERA corrects for inter-gene correlation; a well-known
##     side-effect is conservatism on small designs with strongly
##     co-regulated targets.
##   - ORA hard-thresholds the input at FDR < alpha_fdr_de and
##     loses any chemical whose CTD gene set does not overlap the
##     thresholded list. On small n, this throws away weak true
##     signal -- which is precisely the trade-off ORA is known for
##     and why GSEA was invented.
##
## Reporting Dex's rank under each method is therefore a
## *characterization* of method behavior, not a benchmark of ctdR.
## The script applies BH alpha cutoffs and prints these ranks as-is;
## any method-selection guidance lives in the README / vignette
## "Choosing among the four methods" table.

EXPECTED_HIT <- list(id = "D003907", name = "Dexamethasone")

rank_in <- function(df, id) {
    if (!nrow(df) || !"ChemicalID" %in% colnames(df)) {
        return(data.frame(rank = NA_integer_, total = nrow(df),
            PValue = NA_real_, PValueAdjusted = NA_real_))
    }
    df <- df[order(df$PValue), , drop = FALSE]
    idx <- which(df$ChemicalID == id)
    if (!length(idx)) {
        return(data.frame(rank = NA_integer_, total = nrow(df),
            PValue = NA_real_, PValueAdjusted = NA_real_))
    }
    data.frame(rank = idx[1], total = nrow(df),
        PValue = df$PValue[idx[1]],
        PValueAdjusted = df$PValueAdjusted[idx[1]])
}

hit <- rbind(
    cbind(Method = "ORA",    rank_in(ora_full, EXPECTED_HIT$id)),
    cbind(Method = "GSEA",   rank_in(gsea_full, EXPECTED_HIT$id)),
    cbind(Method = "CAMERA", rank_in(camera_full, EXPECTED_HIT$id))
)
utils::write.table(hit,
    file.path(CONFIG$out_dir, "expected_hit_ranks.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
message("--- Expected-hit check: ", EXPECTED_HIT$name,
    " (", EXPECTED_HIT$id, ") ---")
print(hit, row.names = FALSE)
message("")

## ---- Step G -- plots ---------------------------------------

save_plot <- function(name, plotter) {
    out <- file.path(CONFIG$out_dir, paste0(name, ".png"))
    grDevices::png(out, width = 1600, height = 1200, res = 200)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plotter())
}

## Plot top-15 by raw p-value from the full result so a chart is
## always produced; the alpha-thresholded TSV files (`*_significant`)
## remain the authoritative significance call.
TOP_N_PLOT <- 15L
if (nrow(ora_full)) {
    save_plot("plot_ora",
        function() plot_CTD(ora_full, type = "bar", n = TOP_N_PLOT))
}
if (nrow(gsea_full)) {
    save_plot("plot_gsea",
        function() plot_CTD(gsea_full, type = "bar", n = TOP_N_PLOT))
}
if (nrow(camera_full)) {
    save_plot("plot_camera",
        function() plot_CTD(camera_full, type = "bar", n = TOP_N_PLOT))
}
save_plot("plot_gsva", function() plot_CTD(gsva_scores))
log_step("Step G -- plots written to ", CONFIG$out_dir)

## ---- Step H -- session info --------------------------------

writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(CONFIG$out_dir, "sessionInfo.txt")
)
log_step("Done. Outputs in ", CONFIG$out_dir)
