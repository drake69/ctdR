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

    ## Seed of the random number generator set before GSEA. fgsea
    ## estimates p-values by sampling, so without a seed two runs on the
    ## same input differ: up to 0.9 in a single p-value and 266 against
    ## 271 chemicals at FDR < 0.05 on this dataset. Enrichment scores do
    ## not depend on it.
    gsea_seed = 20261027L,

    out_dir = file.path(getwd(), "example_outputs")
)

dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)

#' Print a progress message prefixed with the time of day.
#'
#' The full pipeline downloads data and runs four enrichment methods on
#' tens of thousands of chemicals; a timestamp on each step shows where
#' the time goes and where a run stopped.
#'
#' @param ... Pieces of the message, pasted together without separator.
log_step <- function(...) {
    message(sprintf("[%s] %s",
        format(Sys.time(), "%H:%M:%S"),
        paste0(..., collapse = "")
    ))
}

## ---- Step A -- a-priori power analysis ---------------------
##
## Question: if a gene truly doubles under Dex (log2FC = 1) and its
## expression varies between samples of the same group with sd = 0.5
## on the log2 scale, how often does a two-sample t-test on that gene
## alone, two-sided at alpha = 0.05, call it significant with 4 DMSO
## and 3 Dex samples? The effect is two standard deviations (1 / 0.5),
## large by any standard; the samples are few.
##
## power.t.test() assumes the two groups have the same size; ours have
## 4 and 3. What drives power is how uncertain the difference between
## the two group means is: sd * sqrt(1/4 + 1/3). Two equal groups of
## 3.43 samples give exactly the same uncertainty (3.43 is the harmonic
## mean of 4 and 3), so we pass n = 3.43.
## The only thing this shortcut gets wrong is the degrees of freedom of
## the t distribution: 4.86 instead of 5. Fewer degrees of freedom make
## the test slightly harder to pass, so the power is slightly
## underestimated (0.554 against an exact 0.559), never overestimated.
## Step D uses limma, not a plain t-test: limma stabilises each gene's
## variance using all the other genes and usually has more power, so
## read 55% as a cautious estimate.

n_dmso_samples <- 4L
n_dex_samples  <- 3L
effective_n_per_group <- 2 / (1 / n_dmso_samples + 1 / n_dex_samples)

power_result <- stats::power.t.test(
    n     = effective_n_per_group,
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
    sprintf("Design : %d DMSO, %d Dex (effective n per group = %.2f, harmonic mean)",
        n_dmso_samples, n_dex_samples, effective_n_per_group),
    sprintf("Target : log2FC = %.2f, residual sd = %.2f",
        CONFIG$power_log2fc_target, CONFIG$power_sd_target),
    sprintf("Alpha  : %.3f (two-sided, per-gene)", CONFIG$alpha_nominal_de),
    sprintf("Power  : %.3f", power_result$power),
    sprintf(""),
    sprintf("Interpretation: ~%.0f%% probability of declaring a true ",
        power_result$power * 100),
    sprintf("log2FC=%.1f gene significant before multiple-testing ",
        CONFIG$power_log2fc_target),
    sprintf("correction. After BH at FDR <= %.2f the realized ",
        CONFIG$alpha_fdr_de),
    sprintf("power per gene is lower; enrichment methods recover ranking "),
    sprintf("signal even when most individual genes do not pass FDR.")
)
writeLines(power_report, file.path(CONFIG$out_dir, "power_report.txt"))
log_step("Step A -- power analysis: power = ",
    sprintf("%.3f", power_result$power))

## ---- Step B -- download + load full GSE311566 --------------

counts_gz_path <- file.path(tempdir(),
    "GSE311566_PBMCs_Female_normalized_counts.txt.gz")
if (!file.exists(counts_gz_path)) {
    log_step("Step B -- downloading ", CONFIG$geo_url)
    utils::download.file(CONFIG$geo_url, counts_gz_path,
        mode = "wb", quiet = TRUE)
}

normalised_counts <- utils::read.table(
    gzfile(counts_gz_path),
    header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
)
normalised_counts$DESCRIPTION <- NULL

## The file also holds the PFAS arms; keep vehicle and Dex only.
dex_vs_dmso_columns <- grep("^(Ctrl|DEX)_F_", colnames(normalised_counts),
    value = TRUE)
stopifnot(length(dex_vs_dmso_columns) >= 6)
normalised_counts <- normalised_counts[, dex_vs_dmso_columns]

treatment_group <- factor(
    ifelse(grepl("^DEX_", dex_vs_dmso_columns), "Dex", "DMSO"),
    levels = c("DMSO", "Dex")
)
samples_per_group <- table(treatment_group)
log_step("Step B -- loaded ", nrow(normalised_counts), " genes x ",
    ncol(normalised_counts), " samples (",
    paste(samples_per_group, names(samples_per_group),
        sep = " ", collapse = ", "),
    ")")

## ---- Step C -- map ENSG -> Entrez --------------------------

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
log_step("Step C -- mapped to Entrez: ", nrow(normalised_counts), " genes")

log2_expression <- log2(as.matrix(normalised_counts) + 1)

## ---- Step D -- differential expression with limma ----------
##
## Data are already library-normalised counts. We fit lmFit + eBayes
## on the log2(x+1) matrix. For raw-count data prefer limma-voom or
## DESeq2; here we follow the GEO supplement's already-normalised
## form to keep the example single-file.

design_matrix <- model.matrix(~ treatment_group)
limma_fit <- limma::lmFit(log2_expression, design_matrix)
limma_fit <- limma::eBayes(limma_fit)
limma_table <- limma::topTable(limma_fit, coef = "treatment_groupDex",
    number = Inf, sort.by = "P")

de_table <- data.frame(
    EntrezID = rownames(limma_table),
    log2FC   = limma_table$logFC,
    pvalue   = limma_table$P.Value,
    padj     = stats::p.adjust(limma_table$P.Value,
        method = CONFIG$p_adjust_method),
    row.names = NULL,
    stringsAsFactors = FALSE
)
utils::write.table(de_table, file.path(CONFIG$out_dir, "de_full.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

n_genes_nominal <- sum(de_table$pvalue < CONFIG$alpha_nominal_de)
n_genes_fdr     <- sum(de_table$padj   < CONFIG$alpha_fdr_de)
log_step(sprintf(
    "Step D -- DE: %d genes pvalue < %.2f, %d genes BH-padj < %.2f",
    n_genes_nominal, CONFIG$alpha_nominal_de,
    n_genes_fdr, CONFIG$alpha_fdr_de
))

## ---- Step E -- verify CTD cache is populated with real data ---

## Ask the package where its cache is and what is in it, rather than
## rebuilding the path here. This check used to call
## rappdirs::user_cache_dir("ctdR"), the location ctdR used before it
## moved to BiocFileCache, so it was validating a directory the package
## no longer writes to: it passed on machines that still had the old
## files lying about and refused to run on clean ones, in both cases
## saying nothing about the data the analysis would actually read.
cached_chemicals <- tryCatch(ctdR::ctd_cache("chemicals"),
    error = function(cache_error) {
        stop(
            conditionMessage(cache_error), "\n",
            "Download CTD_chem_gene_ixns.csv.gz from\n",
            "  https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz\n",
            "gunzip and run:\n",
            "  ctdR::import_CTD(\"path/to/CTD_chem_gene_ixns.csv\")\n",
            "This script intentionally refuses the bundled toy sample.",
            call. = FALSE
        )
    })

## The bundled toy sample has 10 chemicals; the real CTD release has
## tens of thousands. Refuse anything that looks like the toy cache,
## so we never silently report "0 chemicals at FDR < 0.05" coming
## from a 10-chemical universe.
MIN_REAL_CHEMICALS <- 1000L
n_cached_chemicals <- nrow(cached_chemicals)
if (n_cached_chemicals < MIN_REAL_CHEMICALS) {
    stop(
        "The CTD cache holds only ", n_cached_chemicals,
        " chemicals -- this looks like the bundled toy sample.\n",
        "This script requires the real CTD release. Download\n",
        "  https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz\n",
        "gunzip and re-import:\n",
        "  ctdR::import_CTD(\"path/to/CTD_chem_gene_ixns.csv\")\n",
        "before re-running.",
        call. = FALSE
    )
}
ctd_release_info <- tryCatch(ctdR::ctd_provenance(),
    warning = function(provenance_warning) NULL)
log_step("Step E -- CTD cache OK (", n_cached_chemicals, " chemicals",
    if (!is.null(ctd_release_info) &&
        !is.na(ctd_release_info$report_created))
        paste0(", CTD release ", ctd_release_info$report_created) else "",
    ")")

## ---- Step F -- enrichment with all four methods ------------

#' Keep the chemicals whose adjusted p-value passes the threshold.
#'
#' Every method writes both its full table and the significant subset;
#' one helper keeps the significance rule identical across methods.
#'
#' @param enrichment_table Data frame returned by enrichment_CTD(), with
#'   a PValueAdjusted column. Rows with NA adjusted p-value are dropped.
#' @param alpha Numeric in (0, 1), adjusted p-value threshold; a row is
#'   kept when its PValueAdjusted is strictly below it.
#' @return The rows of enrichment_table passing the threshold, possibly
#'   none, with the same columns.
keep_significant <- function(enrichment_table, alpha) {
    passes_threshold <- !is.na(enrichment_table$PValueAdjusted) &
        enrichment_table$PValueAdjusted < alpha
    enrichment_table[passes_threshold, , drop = FALSE]
}

#' Turn list columns into comma-separated strings.
#'
#' fgsea returns the leading edge as a list column, which
#' write.table() cannot write to a TSV.
#'
#' @param enrichment_table Data frame, possibly with list columns.
#' @return The same data frame with every list column replaced by a
#'   character column, elements joined by ",".
flatten_list_columns <- function(enrichment_table) {
    for (column_name in colnames(enrichment_table)) {
        if (is.list(enrichment_table[[column_name]])) {
            enrichment_table[[column_name]] <- vapply(
                enrichment_table[[column_name]],
                function(cell_values) paste(cell_values, collapse = ","),
                character(1))
        }
    }
    enrichment_table
}

#' Write one enrichment table to the output directory as TSV.
#'
#' Keeps file naming and TSV format identical across the four methods.
#'
#' @param enrichment_table Data frame to write; list columns are
#'   flattened first.
#' @param file_stem Character, file name without extension, e.g.
#'   "ora_full"; ".tsv" is appended and the file goes to CONFIG$out_dir.
write_enrichment_table <- function(enrichment_table, file_stem) {
    utils::write.table(flatten_list_columns(enrichment_table),
        file.path(CONFIG$out_dir, paste0(file_stem, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE)
}

## ORA. Choose which column the threshold applies to, then hand over the
## whole table: the genes under the threshold and the background come out
## of the same object, so they cannot disagree.
ora_threshold_column <- "padj"
ora_threshold <- CONFIG$alpha_fdr_de
if (!any(de_table$padj < ora_threshold, na.rm = TRUE)) {
    log_step("Step F -- no genes pass FDR for ORA; ",
        "falling back to nominal p < ", CONFIG$alpha_nominal_de)
    ora_threshold_column <- "pvalue"
    ora_threshold <- CONFIG$alpha_nominal_de
}
# The whole table goes in, with the threshold, rather than a list
# filtered beforehand. The background is then every gene that entered the
# differential test, derived from the same object as the gene list and
# unable to disagree with it. Handing over a pre-filtered list and a
# separate background is the same analysis with one more chance to get it
# wrong: on this dataset the wrong background gives 32 significant
# chemicals against 19, so thirteen of them would be manufactured.
ora_full <- enrichment_CTD(de_table, method = "ORA",
    alpha = ora_threshold, alpha_column = ora_threshold_column,
    pAdjustMethod = CONFIG$p_adjust_method)
ora_significant <- keep_significant(ora_full, CONFIG$alpha_fdr_chemical)
write_enrichment_table(ora_full, "ora_full")
write_enrichment_table(ora_significant, "ora_significant")
log_step("Step F -- ORA: ", nrow(ora_significant),
    " chemicals at FDR < ", CONFIG$alpha_fdr_chemical,
    " (", nrow(ora_full), " total returned)")

## GSEA: full ranked list
set.seed(CONFIG$gsea_seed)
gsea_full <- enrichment_CTD(de_table[, c("EntrezID", "pvalue")],
    method = "GSEA",
    pAdjustMethod = CONFIG$p_adjust_method)
gsea_significant <- keep_significant(gsea_full, CONFIG$alpha_fdr_chemical)
write_enrichment_table(gsea_full, "gsea_full")
write_enrichment_table(gsea_significant, "gsea_significant")
log_step("Step F -- GSEA: ", nrow(gsea_significant),
    " chemicals at FDR < ", CONFIG$alpha_fdr_chemical,
    " (", nrow(gsea_full), " total returned)")

## CAMERA: expression matrix + design + contrast (col 2 = Dex vs DMSO)
camera_full <- enrichment_CTD(log2_expression, method = "CAMERA",
    design = design_matrix, contrast = 2,
    pAdjustMethod = CONFIG$p_adjust_method)
camera_significant <- keep_significant(camera_full,
    CONFIG$alpha_fdr_chemical)
write_enrichment_table(camera_full, "camera_full")
write_enrichment_table(camera_significant, "camera_significant")
log_step("Step F -- CAMERA: ", nrow(camera_significant),
    " chemicals at FDR < ", CONFIG$alpha_fdr_chemical,
    " (", nrow(camera_full), " total returned)")

## GSVA: per-sample scores (no alpha; descriptive output)
gsva_scores <- enrichment_CTD(log2_expression, method = "GSVA")
utils::write.table(
    data.frame(ChemicalID = rownames(gsva_scores), gsva_scores,
        check.names = FALSE),
    file.path(CONFIG$out_dir, "gsva_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
log_step("Step F -- GSVA: ", nrow(gsva_scores), " chemicals x ",
    ncol(gsva_scores), " samples")

## ---- Step F.5 -- top-N summary (printed + saved) -----------

TOP_N <- 10L

#' Take the TOP_N chemicals of one method, ranked by raw p-value.
#'
#' Puts the three tested methods side by side in one table, with the
#' effect size next to the p-value.
#'
#' @param enrichment_table Data frame returned by enrichment_CTD() for
#'   ORA, GSEA or CAMERA. Count, BackgroundRatio and FoldEnrichment are
#'   used when present (ORA) and reported as NA otherwise.
#' @param method_name Character, label written in the Method column.
#' @return Data frame with at most TOP_N rows, or NULL when the method
#'   returned no chemical.
summarise_top_chemicals <- function(enrichment_table, method_name) {
    if (!nrow(enrichment_table)) return(NULL)
    enrichment_table <- enrichment_table[order(enrichment_table$PValue), ,
        drop = FALSE]
    top_rows <- seq_len(min(TOP_N, nrow(enrichment_table)))
    # Effect size travels with significance. Ranking is by p-value, which
    # is the quantity error is controlled on, but a p-value in ORA mixes
    # overlap with set size: a chemical with 16,000 target genes reaches
    # significance at a fold enrichment of 1.35, while one with 2,300
    # needs 2.68. Reporting the rank without the effect hides that.
    input_genes_in_set <- if ("Count" %in% names(enrichment_table))
        enrichment_table$Count[top_rows] else NA_integer_
    chemical_set_size <- if ("BackgroundRatio" %in% names(enrichment_table))
        as.integer(sub("/.*", "", enrichment_table$BackgroundRatio[top_rows]))
        else NA_integer_
    fold_enrichment <- if ("FoldEnrichment" %in% names(enrichment_table))
        enrichment_table$FoldEnrichment[top_rows] else NA_real_
    data.frame(
        Method = method_name,
        Rank = top_rows,
        ChemicalID = enrichment_table$ChemicalID[top_rows],
        ChemicalName = enrichment_table$ChemicalName[top_rows],
        Hits = input_genes_in_set,
        SetSize = chemical_set_size,
        FoldEnrichment = fold_enrichment,
        PValue = enrichment_table$PValue[top_rows],
        PValueAdjusted = enrichment_table$PValueAdjusted[top_rows],
        stringsAsFactors = FALSE
    )
}
top_chemicals <- do.call(rbind, list(
    summarise_top_chemicals(ora_full, "ORA"),
    summarise_top_chemicals(gsea_full, "GSEA"),
    summarise_top_chemicals(camera_full, "CAMERA")
))
utils::write.table(top_chemicals,
    file.path(CONFIG$out_dir, "top_chemicals_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
# Print one block per method, with fixed-width columns, instead of one
# wide data frame. print() on a frame this wide wraps it into three
# detached chunks: every ChemicalID, then every ChemicalName, then every
# p-value, so reading off which chemical ranks third means counting rows
# across three blocks. Chemical names run to seventy characters, so this
# is not an edge case, it is every run.
CHEMICAL_NAME_WIDTH <- 34L

#' Print the top chemicals of one method as a fixed-width block.
#'
#' @param top_chemicals_table Data frame built by
#'   summarise_top_chemicals() for all methods.
#' @param method_name Character, the method whose rows are printed.
#' @param alpha Numeric in (0, 1), adjusted p-value threshold; rows
#'   below it are marked with "*".
#' @return NULL, invisibly; called for the messages it prints.
print_top_chemicals <- function(top_chemicals_table, method_name, alpha) {
    method_rows <- top_chemicals_table[
        top_chemicals_table$Method == method_name, , drop = FALSE]
    message("\n  ", method_name, " -- top ", nrow(method_rows),
        " by raw p-value")
    if (!nrow(method_rows)) {
        message("    (no results)")
        return(invisible(NULL))
    }
    message(sprintf("    %-4s  %-*s  %12s  %5s  %9s  %9s  %s",
        "rank", CHEMICAL_NAME_WIDTH, "chemical", "hits/set", "fold",
        "p", "p.adj", "sig"))
    for (i in seq_len(nrow(method_rows))) {
        chemical_label <- method_rows$ChemicalName[i]
        if (is.na(chemical_label)) chemical_label <- method_rows$ChemicalID[i]
        if (nchar(chemical_label) > CHEMICAL_NAME_WIDTH)
            chemical_label <- paste0(
                substr(chemical_label, 1, CHEMICAL_NAME_WIDTH - 3), "...")
        hits_over_set <- if (is.na(method_rows$Hits[i]) ||
            is.na(method_rows$SetSize[i])) "-"
            else paste0(method_rows$Hits[i], "/", method_rows$SetSize[i])
        fold_label <- if (is.na(method_rows$FoldEnrichment[i])) "-"
            else sprintf("%.2f", method_rows$FoldEnrichment[i])
        message(sprintf("    %-4d  %-*s  %12s  %5s  %9.2e  %9.2e  %s",
            method_rows$Rank[i], CHEMICAL_NAME_WIDTH, chemical_label,
            hits_over_set, fold_label, method_rows$PValue[i],
            method_rows$PValueAdjusted[i],
            if (!is.na(method_rows$PValueAdjusted[i]) &&
                method_rows$PValueAdjusted[i] < alpha) "*" else ""))
    }
    invisible(NULL)
}

message("\n--- Top ", TOP_N, " chemicals per method (by raw p-value) ---")
message("    * marks FDR < ", CONFIG$alpha_fdr_chemical)
message("    hits/set = input genes in the chemical's set / size of that set")
message("    Ranking is by p-value, not by fold enrichment: fold alone has")
message("    no error control, and a two-gene set with both genes hit would")
message("    top it. Read the two together. A large set can reach")
message("    significance at a modest fold, a small one cannot.")
for (method_name in c("ORA", "GSEA", "CAMERA"))
    print_top_chemicals(top_chemicals, method_name, CONFIG$alpha_fdr_chemical)
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

#' Find where one chemical ranks in a method's results.
#'
#' The treatment of the contrast is known, so its rank is a readable
#' check of what each method recovers.
#'
#' @param enrichment_table Data frame returned by enrichment_CTD(),
#'   with ChemicalID, PValue and PValueAdjusted columns.
#' @param chemical_id Character, CTD chemical identifier, e.g. "D003907".
#' @return One-row data frame: rank (1 = smallest raw p-value; NA when
#'   the chemical was not tested), total number of chemicals returned,
#'   and the chemical's PValue and PValueAdjusted (NA when not tested).
rank_of_chemical <- function(enrichment_table, chemical_id) {
    not_tested <- data.frame(rank = NA_integer_,
        total = nrow(enrichment_table),
        PValue = NA_real_, PValueAdjusted = NA_real_)
    if (!nrow(enrichment_table) ||
        !"ChemicalID" %in% colnames(enrichment_table)) {
        return(not_tested)
    }
    enrichment_table <- enrichment_table[order(enrichment_table$PValue), ,
        drop = FALSE]
    matching_rows <- which(enrichment_table$ChemicalID == chemical_id)
    if (!length(matching_rows)) return(not_tested)
    chemical_rank <- matching_rows[1]
    data.frame(rank = chemical_rank, total = nrow(enrichment_table),
        PValue = enrichment_table$PValue[chemical_rank],
        PValueAdjusted = enrichment_table$PValueAdjusted[chemical_rank])
}

expected_hit_ranks <- rbind(
    cbind(Method = "ORA",    rank_of_chemical(ora_full, EXPECTED_HIT$id)),
    cbind(Method = "GSEA",   rank_of_chemical(gsea_full, EXPECTED_HIT$id)),
    cbind(Method = "CAMERA", rank_of_chemical(camera_full, EXPECTED_HIT$id))
)
utils::write.table(expected_hit_ranks,
    file.path(CONFIG$out_dir, "expected_hit_ranks.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
message("--- Expected-hit check: ", EXPECTED_HIT$name,
    " (", EXPECTED_HIT$id, ") ---")
# One line per method, and say in words what a rank of NA means: that the
# chemical never entered the test, which is not the same as having been
# tested and found unremarkable.
for (i in seq_len(nrow(expected_hit_ranks))) {
    hit_rank <- expected_hit_ranks$rank[i]
    if (is.na(hit_rank)) {
        message(sprintf(
            paste0("  %-7s NOT TESTED -- excluded before the test, ",
                "not reported as non-significant (%s chemicals tested)"),
            expected_hit_ranks$Method[i],
            format(expected_hit_ranks$total[i], big.mark = ",")))
    } else {
        message(sprintf("  %-7s rank %s of %s   p = %.2e   p.adj = %.2e%s",
            expected_hit_ranks$Method[i], format(hit_rank, big.mark = ","),
            format(expected_hit_ranks$total[i], big.mark = ","),
            expected_hit_ranks$PValue[i],
            expected_hit_ranks$PValueAdjusted[i],
            if (expected_hit_ranks$PValueAdjusted[i] <
                CONFIG$alpha_fdr_chemical) "   SIGNIFICANT" else ""))
    }
}
message("")

## ---- Step G -- plots ---------------------------------------

#' Render one plot to a PNG in the output directory.
#'
#' The device is opened and closed in one place, so a failing plot
#' cannot leave it open for the next one.
#'
#' @param file_stem Character, file name without extension, e.g.
#'   "plot_ora"; ".png" is appended and the file goes to CONFIG$out_dir.
#' @param build_plot Function with no arguments returning the plot
#'   object to print (1600 x 1200 px at 200 dpi).
save_plot <- function(file_stem, build_plot) {
    png_path <- file.path(CONFIG$out_dir, paste0(file_stem, ".png"))
    grDevices::png(png_path, width = 1600, height = 1200, res = 200)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(build_plot())
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
