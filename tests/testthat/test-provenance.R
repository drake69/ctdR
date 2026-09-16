# The reviewer's point: an analysis is only reproducible if the release of
# CTD behind it can be named. CTD re-releases continuously and does not
# version its filenames, so the date inside the header is the only handle.
# These tests check that the handle reaches the object a user reports on.

.prov_fixture <- function() {
    cols <- c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID", "Interaction",
        "InteractionActions", "PubMedIDs")
    rows <- c(
        "Acetaminophen,D000082,,TNF,7124,,Homo sapiens,9606,x,increases^expression,1",
        "Acetaminophen,D000082,,IL6,3569,,Homo sapiens,9606,x,increases^expression,2",
        "Benzene,D001554,,TP53,7157,,Homo sapiens,9606,x,affects^binding,3"
    )
    tmp <- tempfile(fileext = ".csv")
    writeLines(c(
        "# preamble",
        "# Report created: Thu May 28 13:30:51 EDT 2026",
        "#", "# Fields:", paste0("# ", paste(cols, collapse = ",")), "#",
        rows
    ), tmp)
    tmp
}

.with_temp_cache <- function(code) {
    tmp_cache <- file.path(tempdir(), paste0("ctdR_prov_", basename(tempfile())))
    options(ctdR.cache = tmp_cache)
    on.exit({
        .restore_ctd_cache()
        unlink(tmp_cache, recursive = TRUE)
    }, add = TRUE)
    force(code)
}

test_that("import_CTD records the CTD release date from the file header", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    f <- .prov_fixture()
    on.exit(unlink(f))

    .with_temp_cache({
        suppressMessages(suppressWarnings(import_CTD(f)))
        genes <- data.frame(EntrezID = c("7124", "3569"), pvalue = c(.01, .02))
        res <- suppressMessages(suppressWarnings(
            enrichment_CTD(genes, method = "ORA")))
        prov <- ctd_provenance(res)

        expect_s3_class(prov, "ctd_provenance")
        expect_identical(prov$report_created, "Thu May 28 13:30:51 EDT 2026")
        expect_identical(prov$source, f)
        expect_identical(prov$ctdR_version,
            as.character(utils::packageVersion("ctdR")))
        expect_true(prov$n_chemicals >= 1L)
        expect_true(prov$n_interactions >= 1L)
        # The import timestamp is when the data was read, not when CTD
        # built it: both are needed, and they are different questions.
        expect_false(identical(prov$accessed, prov$report_created))
    })
})

test_that("provenance reaches every method, in the right container", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("GSVA")

    f <- .prov_fixture()
    on.exit(unlink(f))

    .with_temp_cache({
        suppressMessages(suppressWarnings(import_CTD(f)))
        genes <- data.frame(EntrezID = c("7124", "3569"), pvalue = c(.01, .02))

        ora <- suppressMessages(suppressWarnings(
            enrichment_CTD(genes, method = "ORA")))
        # A data frame has no metadata slot, so the record rides on an
        # attribute; the accessor hides which is which.
        expect_false(is.null(attr(ora, "ctd_provenance")))
        expect_s3_class(ctd_provenance(ora), "ctd_provenance")

        se <- readRDS(system.file("extdata", "GSE311566_subset.rds",
            package = "ctdR"))
        gsva <- suppressMessages(suppressWarnings(
            enrichment_CTD(se, method = "GSVA")))
        # A SummarizedExperiment does have one, and that is where it goes.
        expect_true(.is_se(gsva))
        expect_false(is.null(S4Vectors::metadata(gsva)$ctd_provenance))
        expect_s3_class(ctd_provenance(gsva), "ctd_provenance")
    })
})

test_that("ctd_provenance() with no argument reads the cache", {
    # Asking which release is about to be analysed comes before any
    # analysis exists to ask it of.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")

    f <- .prov_fixture()
    on.exit(unlink(f))

    .with_temp_cache({
        expect_warning(empty <- ctd_provenance(), "No CTD provenance in the cache")
        expect_null(empty)

        suppressMessages(suppressWarnings(import_CTD(f)))
        prov <- ctd_provenance()
        expect_s3_class(prov, "ctd_provenance")
        expect_identical(prov$report_created, "Thu May 28 13:30:51 EDT 2026")
    })
})

test_that("provenance survives subsetting and ordering", {
    prov <- .ctd_provenance_new("Thu May 28 13:30:51 EDT 2026", "f.csv", 3L, 9L)
    df <- .attach_provenance(
        data.frame(ChemicalID = c("A", "B", "C"), PValue = c(.3, .1, .2)), prov)

    expect_s3_class(ctd_provenance(df[1:2, ]), "ctd_provenance")
    expect_s3_class(ctd_provenance(head(df, 2)), "ctd_provenance")
    expect_s3_class(ctd_provenance(df[order(df$PValue), ]), "ctd_provenance")
    expect_s3_class(ctd_provenance(rbind(df, df)), "ctd_provenance")
})

test_that("ctd_provenance warns rather than returning a silent NULL", {
    # merge() drops attributes. Reporting a blank provenance would be worse
    # than saying so, because the user is documenting an analysis with it.
    prov <- .ctd_provenance_new("Thu May 28 13:30:51 EDT 2026", "f.csv", 3L, 9L)
    df <- .attach_provenance(
        data.frame(ChemicalID = c("A", "B"), PValue = c(.1, .2)), prov)

    merged <- merge(df, df, by = "ChemicalID")
    expect_warning(out <- ctd_provenance(merged), "No CTD provenance")
    expect_null(out)
})

test_that("a file without a Report created line still yields a record", {
    # The date is CTD's to provide. When it is absent the rest of the
    # record is still worth having, and the gap is stated rather than
    # papered over with the import date.
    prov <- .ctd_provenance_new(NA_character_, "f.csv", 3L, 9L)

    expect_true(is.na(prov$report_created))
    expect_output(print(prov), "not stated in the file")
})

test_that("the printed record names the release and the import separately", {
    prov <- .ctd_provenance_new("Thu May 28 13:30:51 EDT 2026", "f.csv", 3L, 9L)

    expect_output(print(prov), "CTD provenance")
    expect_output(print(prov), "Report created: Thu May 28 13:30:51 EDT 2026")
    expect_output(print(prov), "Imported:")
    expect_output(print(prov), "3 chemicals, 9 chemical-gene pairs")
})

test_that("a cache written before provenance existed still works", {
    # Users who imported with an earlier version keep a cache with no
    # provenance entry. That must degrade to "no record", not to an error.
    expect_identical(.attach_provenance(data.frame(a = 1), NULL),
        data.frame(a = 1))
})
