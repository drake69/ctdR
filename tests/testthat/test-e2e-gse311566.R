## End-to-end biological-expectation test on the bundled GSE311566 subset.
##
## Mirrors the vignette pipeline:
##   1. import the toy CTD sample
##   2. load inst/extdata/GSE311566_subset.rds (Dex vs DMSO, female PBMCs)
##   3. base-R differential expression (per-gene t.test + BH)
##   4. run GSEA and CAMERA via enrichment_CTD()
##
## Guards against regressions in ID mapping, schema, or sorting that the
## vignette/example chunks would only catch as "still runs", not as
## "still returns the biologically expected hit".
##
## Expected: Dexamethasone (CTD ID D003907) surfaces near the top in both
## methods. Thresholds are deliberately loose to absorb fgsea's tie-breaking
## non-determinism and t.test noise on a 4-vs-3 sample design.

test_that("e2e GSE311566 Dex-vs-DMSO recovers Dexamethasone via GSEA + CAMERA", {
    skip_on_cran()
    skip_if_not_installed("fgsea")
    skip_if_not_installed("limma")
    skip_if_not_installed("org.Hs.eg.db")

    subset_path <- system.file(
        "extdata", "GSE311566_subset.rds", package = "ctdR"
    )
    ctd_sample <- system.file(
        "extdata", "CTD_chem_gene_ixns_sample.csv", package = "ctdR"
    )
    skip_if(subset_path == "" || ctd_sample == "",
        "bundled GSE311566 / CTD sample not found in installed package")

    tmp_cache <- file.path(tempdir(), "ctdR_e2e_test")
    dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)
    options(ctdR.cache = tmp_cache)
    on.exit({
        options(ctdR.cache = NULL)
        unlink(tmp_cache, recursive = TRUE)
    })

    suppressMessages(import_CTD(ctd_sample))

    se <- readRDS(subset_path)
    ## The bundled object is a SummarizedExperiment: assay and sample
    ## annotation travel together, so a subset or a reorder cannot
    ## desynchronise the group labels from the columns they describe.
    expect_s4_class(se, "SummarizedExperiment")
    expr <- SummarizedExperiment::assay(se)
    grp  <- se$group
    expect_identical(levels(grp), c("DMSO", "Dex"))
    expect_equal(ncol(se), 7L)
    expect_identical(colnames(se), colnames(expr))

    de <- t(apply(expr, 1, function(y) {
        tt <- stats::t.test(y ~ grp)
        c(pvalue = tt$p.value)
    }))
    de <- data.frame(
        EntrezID = rownames(expr),
        pvalue = as.numeric(de),
        stringsAsFactors = FALSE
    )

    ## --- GSEA: Dex must be among the top 3 by raw p-value -----------------
    suppressWarnings({
        gsea_res <- enrichment_CTD(de, method = "GSEA")
    })
    expect_s3_class(gsea_res, "data.frame")
    expect_true(all(c("ChemicalID", "ChemicalName", "Method",
                      "PValue", "PValueAdjusted") %in% colnames(gsea_res)))
    expect_equal(unique(gsea_res$Method), "GSEA")

    gsea_by_p <- gsea_res[order(gsea_res$PValue), ]
    expect_true("D003907" %in% head(gsea_by_p$ChemicalID, 3),
        info = "Dexamethasone (D003907) should rank in the top 3 by GSEA p-value")
    # No absolute PValue threshold: fgsea tie-breaking is non-deterministic
    # and the 4-vs-3 t.test stats vector has many ties on this small subset.
    # The rank-in-top-3 assertion above is the meaningful biological guard.

    ## --- CAMERA: Dex must be among the top 6 (all 10 ranked) --------------
    ## Fed the SummarizedExperiment directly, not its assay.
    design <- stats::model.matrix(~ grp)
    cam_res <- enrichment_CTD(se, method = "CAMERA",
        design = design, contrast = 2)
    expect_s3_class(cam_res, "data.frame")
    expect_true("Direction" %in% colnames(cam_res))
    cam_by_p <- cam_res[order(cam_res$PValue), ]
    expect_true("D003907" %in% head(cam_by_p$ChemicalID, 6),
        info = "Dexamethasone (D003907) should rank in the top 6 by CAMERA p-value")

    ## --- GSVA: SummarizedExperiment in, SummarizedExperiment out ----------
    skip_if_not_installed("GSVA")
    suppressMessages({
        gsva_res <- enrichment_CTD(se, method = "GSVA")
    })
    expect_s4_class(gsva_res, "SummarizedExperiment")
    expect_equal(ncol(gsva_res), 7L)
    expect_true("D003907" %in% rownames(gsva_res))
    ## Sample annotation survives the round trip, which is the point of
    ## accepting the container rather than the bare matrix.
    expect_identical(gsva_res$group, grp)
    expect_identical(colnames(gsva_res), colnames(se))
})
