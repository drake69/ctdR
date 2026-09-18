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
#' Earlier versions delegated this step to the \pkg{clusterProfiler}
#' package. That pulled in 59 packages, and a visualization layer this
#' package never used, to reach a single call. The test itself is one
#' line of \code{stats}, so it is computed here and that dependency is
#' gone. The p-values are
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
#' @param universe Vector of background gene identifiers, or \code{NULL}
#'   (default) to use every gene in \code{ChemicalName_GeneSymbols}. Set
#'   it to \code{rownames(expr)} or to the full tested gene list to
#'   restrict the background to measured genes only. Coerced with
#'   \code{as.character()}, so an integer column of Entrez IDs works as
#'   it does for \code{gene_symbols}.
#'
#'   The default is a fallback, not a recommendation, and the function
#'   says so when it uses it. A hypergeometric test asks how many of
#'   \eqn{n} genes drawn from \eqn{N} would land in a set. Genes that
#'   your experiment could never have detected still sit in \eqn{N},
#'   filling the urn with balls that cannot be drawn, so the draw looks
#'   more selective than it was and the p-value comes out too small. The
#'   error is anti-conservative: it manufactures significance.
#'
#'   The right universe is every gene that entered your test, not every
#'   gene you sequenced and not only the significant ones: a gene
#'   filtered out for low expression could not have been selected, so it
#'   does not belong in the urn either. On the RNA-seq example bundled
#'   with this package the difference is 32 significant chemicals
#'   against 19.
#' @param minGSSize Integer. Minimum gene set size after intersection
#'   with the background (default 2). The default is chosen for CTD,
#'   where the median chemical has 4 target genes: a one-gene set is
#'   degenerate, since its p-value equals the ratio of input genes to
#'   background whichever gene it contains, so it measures membership
#'   rather than enrichment.
#' @param maxGSSize Maximum gene set size after intersection with the
#'   background. The default is \code{Inf}, that is no upper limit, and
#'   that too is a choice made for CTD rather than inherited.
#'
#'   The reasoning mirrors the one for \code{minGSSize}, and reaches the
#'   opposite conclusion. A set of \eqn{M} genes cannot, even when every
#'   one of the \eqn{m} input genes falls inside it, produce a p-value
#'   below \eqn{C(M,m)/C(N,m) \approx (M/N)^m}. That floor rises with
#'   \eqn{M}, so a large enough set becomes untestable. Measured on CTD,
#'   with 28,571 genes in the background and an input list of 169, the
#'   largest still-testable set is about 26,600 genes, 93\% of the
#'   universe. The largest chemical in CTD has 16,536. No chemical is
#'   untestable from above, so an upper cut removes sets that could have
#'   been declared significant.
#'
#'   What it removes is not marginal. A cut at 500 excludes 265
#'   chemicals, among them benzo(a)pyrene, valproic acid, sodium
#'   arsenite, bisphenol A, aflatoxin B1 and particulate matter. Their
#'   sets are large because the literature on them is large: in CTD size
#'   tracks how well studied a chemical is, where in GO a large term is
#'   one that has stopped meaning anything. Keeping them costs 3\% more
#'   tests, 8,235 against 7,970.
#'
#'   Set it to a finite value if you have a reason of your own; ctdR
#'   does not impose one.
#'
#'   One consequence to be aware of: without a cap the significant
#'   results skew towards large sets. This is not the fold enrichment
#'   talking, which moves the other way, since \eqn{M} sits in its
#'   denominator. It is that a p-value measures how unlikely an excess
#'   is, and the excess is counted in genes. With 154 input genes from a
#'   background of 27,444, a fold of 1.5 means 0.4 genes above
#'   expectation for a set of 100 (p = 0.43) and 45 genes above for a
#'   set of 16,000 (p = 1.4e-15). Sort by \code{p.adjust} to rank by
#'   evidence, and read \code{foldEnrichment} beside it to see how sharp
#'   the association is; neither answers the question alone.
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
    minGSSize = 2, maxGSSize = Inf) {
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
        # Coerce rather than reject. Gene identifiers arrive as integers
        # often enough (a table round-tripped through read.delim, say)
        # that rejecting them would break the most ordinary use of this
        # argument, restricting the background to measured genes, while
        # the same column passed as the input list is already coerced.
        if (!is.atomic(universe))
            stop("'universe' must be a vector of gene identifiers, not a ",
                class(universe)[1], ".", call. = FALSE)
        universe <- as.character(universe)
        background <- intersect(background, universe)
    }
    N <- length(background)
    if (is.null(universe))
        message(sprintf(
            paste0("background: all %d genes in the CTD sets, because no ",
                "'universe' was given. If your experiment could only ",
                "detect some of them, pass those as 'universe': a ",
                "background wider than what was measurable makes ",
                "p-values too small."),
            N))

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
            "gene set size filter [%s, %s]: %d of %d chemicals %s",
            format(minGSSize), format(maxGSSize),
            sum(!keep), length(keep),
            sprintf("not tested (%d below, %d above); %d tested.",
                sum(sizes < minGSSize), sum(sizes > maxGSSize), sum(keep))))
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
