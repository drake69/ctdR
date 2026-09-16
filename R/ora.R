#' @title Over-Representation Analysis (ORA)
#'
#' @description
#' Internal engine that performs Over-Representation Analysis with a
#' hypergeometric test computed directly on \code{\link[stats]{phyper}}.
#' Returns a result table with one row per tested chemical.
#'
#' Column renaming, multiple-testing correction, chemical-name join,
#' canonical ordering and sort are applied downstream by
#' \code{\link{.format_enrichment_result}} so the engine stays close
#' to the vocabulary of the test itself.
#'
#' @details
#' Earlier versions delegated this step to
#' \code{clusterProfiler::enricher()}. That pulled in 59 packages, and a
#' visualization layer this package never used, to reach a single call.
#' The test itself is one line of \code{stats}, so it is computed here
#' and \code{clusterProfiler} is no longer a dependency. The p-values are
#' unchanged: the two implementations were compared over 24
#' configurations of input list, minimum set size and background
#' universe, and the largest absolute difference was exactly 0.
#'
#' The background is every gene appearing in \code{ChemicalName_GeneSymbols},
#' optionally narrowed by \code{universe}. Gene sets are intersected with
#' that background \emph{before} the size filter, so \code{minGSSize} and
#' \code{maxGSSize} always refer to the set as actually tested rather than
#' to its nominal size.
#'
#' @param ChemicalName_GeneSymbols A data frame with two columns
#'   (\code{term}, \code{gene}) mapping CTD chemical IDs to HGNC gene
#'   symbols. The first two columns are used, whatever their names.
#' @param gene_symbols Character vector of HGNC gene symbols to test
#'   for enrichment. Duplicates and \code{NA} are removed.
#' @param pAdjustMethod Character. Method for multiple testing
#'   correction (default \code{"BH"}). Passed to
#'   \code{\link[stats]{p.adjust}}.
#' @param universe Character vector of background gene symbols, or
#'   \code{NULL} (default) to use every gene in
#'   \code{ChemicalName_GeneSymbols}. Set it to \code{rownames(expr)} or
#'   to the full tested gene list to restrict the background to measured
#'   genes only.
#' @param minGSSize Integer. Minimum gene set size after intersection
#'   with the background (default 2). The default is chosen for CTD,
#'   where the median chemical has 4 target genes: a one-gene set is
#'   degenerate, since its p-value equals the ratio of input genes to
#'   background whichever gene it contains, so it measures membership
#'   rather than enrichment.
#' @param maxGSSize Integer. Maximum gene set size after intersection
#'   with the background (default 500).
#'
#' @return A data frame with columns \code{ChemicalID}, \code{GeneRatio},
#'   \code{BgRatio}, \code{pvalue}, \code{p.adjust}, \code{geneID},
#'   \code{Count} and \code{foldEnrichment}, sorted by \code{pvalue}
#'   ascending. Returns an empty data frame with the same structure when
#'   no gene set can be tested. Emits a message reporting how many
#'   chemicals the size filter left untested, since those are absent from
#'   the result rather than present with a large p-value.
#'
#' @keywords internal
ora <- function(ChemicalName_GeneSymbols, gene_symbols,
    pAdjustMethod = "BH", universe = NULL,
    minGSSize = 2, maxGSSize = 500) {
    empty <- data.frame(
        ChemicalID = character(), GeneRatio = character(),
        BgRatio = character(), pvalue = numeric(),
        p.adjust = numeric(), geneID = character(),
        Count = integer(), foldEnrichment = numeric(),
        stringsAsFactors = FALSE
    )

    t2g <- as.data.frame(ChemicalName_GeneSymbols)
    if (ncol(t2g) < 2L)
        stop("'ChemicalName_GeneSymbols' needs a term and a gene column.",
            call. = FALSE)
    term <- as.character(t2g[[1L]])
    gene <- as.character(t2g[[2L]])

    background <- unique(gene[!is.na(gene)])
    if (!is.null(universe)) {
        if (!is.character(universe))
            stop("'universe' must be a character vector.", call. = FALSE)
        background <- intersect(background, universe)
    }
    N <- length(background)

    hits <- unique(as.character(gene_symbols))
    hits <- intersect(hits[!is.na(hits)], background)
    n <- length(hits)
    if (N == 0L || n == 0L) return(empty)

    gene_sets <- split(gene, term)
    gene_sets <- lapply(gene_sets, function(g) intersect(unique(g), background))
    sizes <- lengths(gene_sets)
    keep <- sizes >= minGSSize & sizes <= maxGSSize
    # Report what the filter removed. A chemical dropped here is not a
    # non-significant result, it was never tested, and silently missing
    # rows are indistinguishable from rows that came back empty.
    if (any(!keep))
        message(sprintf(
            paste0("gene set size filter [%d, %d]: %d of %d chemicals ",
                "not tested (%d below, %d above); %d tested."),
            as.integer(minGSSize), as.integer(maxGSSize),
            sum(!keep), length(keep),
            sum(sizes < minGSSize), sum(sizes > maxGSSize), sum(keep)))
    gene_sets <- gene_sets[keep]
    if (!length(gene_sets)) return(empty)

    M <- unname(lengths(gene_sets))
    # Keep the input order inside geneID so the column is reproducible.
    overlap <- lapply(gene_sets, function(g) hits[hits %in% g])
    k <- unname(lengths(overlap))

    # P(X >= k) for X hypergeometric: k or more of the n input genes
    # falling in a set of M, drawn from a background of N.
    pvalue <- stats::phyper(k - 1L, M, N - M, n, lower.tail = FALSE)

    res <- data.frame(
        ChemicalID = names(gene_sets),
        GeneRatio = paste0(k, "/", n),
        BgRatio = paste0(M, "/", N),
        pvalue = pvalue,
        p.adjust = stats::p.adjust(pvalue, method = pAdjustMethod),
        geneID = unname(vapply(overlap, paste, character(1), collapse = "/")),
        Count = as.integer(k),
        foldEnrichment = (k / n) / (M / N),
        stringsAsFactors = FALSE
    )
    res <- res[order(res$pvalue), , drop = FALSE]
    rownames(res) <- NULL
    res
}
