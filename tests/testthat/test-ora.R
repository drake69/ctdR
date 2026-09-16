# The reference for these tests is the hypergeometric distribution
# itself, not another package: ora() computes P(X >= k) directly, so the
# expected values below are stats::phyper() calls and closed-form
# combinatorics rather than a second implementation's output.
#
# Note throughout: the background is the union of the gene sets, not a
# pool supplied separately, so each fixture carries a FILLER set whose
# only job is to pin the background size at a readable number.

test_that("ora matches the hypergeometric tail probability", {
    target <- paste0("TARGET", 1:20)
    bg <- paste0("GENE", 1:300)
    term2gene <- rbind(
        data.frame(term = "CHEM1", gene = target),
        data.frame(term = "CHEM2", gene = bg[1:50]),
        data.frame(term = "CHEM3", gene = bg[201:300]),
        data.frame(term = "FILLER", gene = bg)
    )
    # Background N = 320 (20 targets + 300 others); the input list is the
    # 20 targets plus 10 background genes, so n = 30.
    gene_list <- c(target, bg[201:210])

    result <- ora(term2gene, gene_list)

    expect_true(is.data.frame(result))
    expect_identical(result$ChemicalID[1], "CHEM1")

    # CHEM1: perfect overlap, k = 20 of a set of M = 20.
    expect_equal(
        result$pvalue[result$ChemicalID == "CHEM1"],
        stats::phyper(20 - 1, 20, 320 - 20, 30, lower.tail = FALSE)
    )
    # CHEM3: k = 10 of M = 100, below what chance alone predicts.
    expect_equal(
        result$pvalue[result$ChemicalID == "CHEM3"],
        stats::phyper(10 - 1, 100, 320 - 100, 30, lower.tail = FALSE)
    )
    # CHEM2: no overlap at all, so P(X >= 0) is exactly 1.
    expect_equal(result$pvalue[result$ChemicalID == "CHEM2"], 1)
})

test_that("ora reports ratios, counts and fold enrichment consistently", {
    target <- paste0("TARGET", 1:20)
    bg <- paste0("GENE", 1:300)
    term2gene <- rbind(
        data.frame(term = "CHEM1", gene = target),
        data.frame(term = "FILLER", gene = bg)
    )
    result <- ora(term2gene, c(target, bg[201:210]))
    chem1 <- result[result$ChemicalID == "CHEM1", ]

    expect_identical(chem1$GeneRatio, "20/30")
    expect_identical(chem1$BgRatio, "20/320")
    expect_identical(chem1$Count, 20L)
    # FoldEnrichment is (k/n) / (M/N), here (20/30) / (20/320).
    expect_equal(chem1$foldEnrichment, (20 / 30) / (20 / 320))
    expect_identical(
        sort(strsplit(chem1$geneID, "/", fixed = TRUE)[[1]]),
        sort(target)
    )
    # Results are sorted by p-value, not by chemical ID.
    expect_false(is.unsorted(result$pvalue))
})

test_that("no upper limit applies unless one is asked for", {
    # A set of 600 genes is testable and must be tested by default. The
    # measured ceiling on CTD is around 26,600 genes, well above the
    # largest chemical in the database.
    bg <- paste0("GENE", 1:1000)
    term2gene <- rbind(
        data.frame(term = "BIG", gene = bg[1:600]),
        data.frame(term = "SMALL", gene = bg[1:5]),
        data.frame(term = "FILLER", gene = bg)
    )
    gene_list <- bg[1:20]

    expect_true("BIG" %in% ora(term2gene, gene_list)$ChemicalID)
    expect_false("BIG" %in% ora(term2gene, gene_list,
        maxGSSize = 500)$ChemicalID)
})

test_that("a one-gene set has the same p-value whatever gene it holds", {
    # This is why minGSSize defaults to 2: for M = 1 the hypergeometric
    # p-value collapses to n/N, so the test reports membership rather
    # than enrichment and carries no information about the chemical.
    bg <- paste0("GENE", 1:100)
    term2gene <- rbind(
        data.frame(term = "ONE_A", gene = bg[1]),
        data.frame(term = "ONE_B", gene = bg[2]),
        data.frame(term = "FILLER", gene = bg)
    )
    gene_list <- bg[1:10]

    result <- ora(term2gene, gene_list, minGSSize = 1)
    singles <- result[result$ChemicalID %in% c("ONE_A", "ONE_B"), ]

    expect_identical(nrow(singles), 2L)
    expect_equal(singles$pvalue[1], singles$pvalue[2])
    expect_equal(singles$pvalue[1], 10 / 100)
})

test_that("the size thresholds are chosen for CTD, not inherited", {
    # Both defaults were once taken from a general-purpose tool. The lower
    # one is now 2 because CTD sets are small; the upper one is gone
    # because no CTD set is large enough to be untestable.
    expect_identical(formals(ora)$minGSSize, 2)
    expect_identical(formals(ora)$maxGSSize, Inf)

    bg <- paste0("GENE", 1:100)
    term2gene <- rbind(
        data.frame(term = "SINGLETON", gene = bg[1]),
        data.frame(term = "PAIR", gene = bg[2:3]),
        data.frame(term = "FILLER", gene = bg)
    )
    result <- ora(term2gene, bg[1:10])

    expect_false("SINGLETON" %in% result$ChemicalID)
    expect_true("PAIR" %in% result$ChemicalID)
})

test_that("ora reports the chemicals the size filter left untested", {
    # A chemical dropped by the filter is absent from the result, not
    # present with a large p-value. Without the message the two cases
    # are indistinguishable to the caller.
    bg <- paste0("GENE", 1:100)
    term2gene <- rbind(
        data.frame(term = "SINGLETON_A", gene = bg[1]),
        data.frame(term = "SINGLETON_B", gene = bg[2]),
        data.frame(term = "PAIR", gene = bg[3:4]),
        data.frame(term = "FILLER", gene = bg)
    )

    expect_message(
        ora(term2gene, bg[1:10]),
        "2 of 4 chemicals not tested \\(2 below, 0 above\\); 2 tested"
    )
    # Nothing filtered means nothing to report.
    expect_no_message(ora(term2gene, bg[1:10], minGSSize = 1))
})

test_that("the size filter applies after intersection with the background", {
    # SPARSE nominally has 6 genes but only 2 survive the universe, so a
    # minGSSize of 3 must drop it: the filter describes the set as
    # tested, not as declared.
    term2gene <- rbind(
        data.frame(term = "SPARSE", gene = c("A", "B", "C", "D", "E", "F")),
        data.frame(term = "FILLER", gene = paste0("F", 1:20))
    )
    universe <- c("A", "B", paste0("F", 1:20))

    kept <- ora(term2gene, c("A", "B"), universe = universe, minGSSize = 2)
    dropped <- ora(term2gene, c("A", "B"), universe = universe, minGSSize = 3)

    expect_true("SPARSE" %in% kept$ChemicalID)
    # Background is the intersection of the sets with the universe: 22.
    expect_identical(kept$BgRatio[kept$ChemicalID == "SPARSE"], "2/22")
    expect_false("SPARSE" %in% dropped$ChemicalID)
})

test_that("universe narrows the background and is validated", {
    term2gene <- rbind(
        data.frame(term = "CHEM1", gene = paste0("G", 1:10)),
        data.frame(term = "CHEM2", gene = paste0("G", 11:20))
    )
    # Without a universe the background is every gene in term2gene (20).
    full <- ora(term2gene, paste0("G", 1:5))
    expect_identical(full$BgRatio[full$ChemicalID == "CHEM1"], "10/20")

    narrowed <- ora(term2gene, paste0("G", 1:5), universe = paste0("G", 1:12))
    expect_identical(narrowed$BgRatio[narrowed$ChemicalID == "CHEM1"], "10/12")

    # An integer column of Entrez IDs, which is what a DE table read back
    # with read.delim() gives, must work as a universe: the same column
    # is already coerced when it arrives as the input gene list.
    int_universe <- ora(
        data.frame(term = "CHEM1", gene = as.character(1:10)),
        c("1", "2"), universe = 1:12)
    expect_identical(int_universe$BgRatio, "10/10")

    factor_universe <- ora(term2gene, paste0("G", 1:5),
        universe = factor(paste0("G", 1:12)))
    expect_identical(
        factor_universe$BgRatio[factor_universe$ChemicalID == "CHEM1"],
        "10/12")

    expect_error(ora(term2gene, "G1", universe = list("G1")),
        "vector of gene identifiers")
})

test_that("ora returns the empty schema when nothing can be tested", {
    term2gene <- data.frame(term = rep("CHEM1", 2), gene = c("TP53", "TNF"))
    cols <- c("ChemicalID", "GeneRatio", "BgRatio", "pvalue",
        "p.adjust", "geneID", "Count", "foldEnrichment")

    # No input gene is in the background.
    none <- ora(term2gene, c("NONEXISTENT_1", "NONEXISTENT_2"))
    expect_identical(nrow(none), 0L)
    expect_identical(colnames(none), cols)

    # No gene set survives the size filter.
    too_big <- ora(term2gene, "TP53", minGSSize = 50)
    expect_identical(nrow(too_big), 0L)
    expect_identical(colnames(too_big), cols)
})

test_that("ora drops duplicate and NA input genes before testing", {
    term2gene <- rbind(
        data.frame(term = "CHEM1", gene = paste0("G", 1:10)),
        data.frame(term = "FILLER", gene = paste0("B", 1:90))
    )
    plain <- ora(term2gene, c("G1", "G2", "B1"))
    noisy <- ora(term2gene, c("G1", "G1", "G2", NA, "G2", "B1"))

    expect_equal(noisy$pvalue, plain$pvalue)
    # n counts distinct input genes present in the background, so 3.
    expect_identical(noisy$GeneRatio[noisy$ChemicalID == "CHEM1"], "2/3")
})

test_that("pAdjustMethod reaches stats::p.adjust", {
    bg <- paste0("GENE", 1:200)
    term2gene <- do.call(rbind, lapply(1:10, function(i)
        data.frame(term = paste0("CHEM", i),
            gene = bg[((i - 1) * 20 + 1):(i * 20)])))
    gene_list <- bg[1:25]

    bh <- ora(term2gene, gene_list, pAdjustMethod = "BH")
    bonf <- ora(term2gene, gene_list, pAdjustMethod = "bonferroni")

    expect_equal(bh$p.adjust, stats::p.adjust(bh$pvalue, "BH"))
    expect_equal(bonf$p.adjust, stats::p.adjust(bonf$pvalue, "bonferroni"))
    expect_true(all(bonf$p.adjust >= bh$p.adjust))
})
