# A CTD download has no header row: the field names live inside the
# commented preamble, on the line after "# Fields:". These fixtures mirror
# that shape rather than the shape the old hard-coded skip expected, so
# what the tests exercise is what users actually hand the package.
.write_fake_ctd <- function(col_names, data_rows = character(0),
    n_preamble = 12L, report_created = "Thu May 28 13:30:51 EDT 2026",
    fields_marker = TRUE) {
    preamble <- paste0("# line ", seq_len(n_preamble))
    if (!is.na(report_created))
        preamble <- c(preamble, paste("# Report created:", report_created))
    header <- c(preamble, "#",
        if (fields_marker) "# Fields:",
        paste0("# ", paste(col_names, collapse = ",")), "#")
    tmp <- tempfile(fileext = ".csv")
    writeLines(c(header, data_rows), tmp)
    tmp
}

.ctd_cols <- function() {
    c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID", "Interaction",
        "InteractionActions", "PubMedIDs")
}

test_that("import_CTD errors on non-existent file", {
    expect_error(
        import_CTD("/no/such/file.csv"),
        "File not found"
    )
})

test_that("import_CTD errors with helpful download message", {
    expect_error(
        import_CTD("/no/such/file.csv"),
        "ctdbase.org"
    )
})

test_that("import_CTD errors when the header names no known field", {
    # A file whose commented header lists none of the CTD field names is
    # rejected before any record is read, and the message states the rule
    # that was applied rather than a line number, which would be useless
    # against a file whose shape has changed.
    tmp_file <- .write_fake_ctd(c("WrongCol1", "WrongCol2", "WrongCol3"),
        "val1,val2,val3")
    on.exit(unlink(tmp_file))

    expect_error(import_CTD(tmp_file), "Could not find the column-name line")
    expect_error(import_CTD(tmp_file), "at least 3 of these names")
    expect_error(import_CTD(tmp_file), "ChemicalID")
})

test_that("import_CTD errors on a CTD-shaped file missing required columns", {
    # Enough known names to be recognised as the header, not enough to be
    # a usable chemical-gene file.
    tmp_file <- .write_fake_ctd(
        c("ChemicalName", "ChemicalID", "GeneSymbol", "SomethingElse"),
        c("a,b,c,d", "e,f,g,h"))
    on.exit(unlink(tmp_file))

    expect_error(import_CTD(tmp_file), "Not a valid CTD file")
    expect_error(import_CTD(tmp_file), "Missing columns")
})

test_that("the last line naming the fields wins over an earlier mention", {
    # The preamble of a real CTD file is prose, and prose can name columns.
    # The line that counts is the one just before the data.
    cols <- .ctd_cols()
    tmp <- tempfile(fileext = ".csv")
    writeLines(c(
        "# This file lists ChemicalName, ChemicalID and GeneSymbol per row.",
        "# Report created: Thu May 28 13:30:51 EDT 2026",
        "# Fields:",
        paste0("# ", paste(cols, collapse = ",")),
        "#",
        "Acetaminophen,D000082,,TNF,7124,,Homo sapiens,9606,x,y,1"
    ), tmp)
    on.exit(unlink(tmp))

    hdr <- ctdR:::.parse_ctd_header(tmp)
    expect_identical(hdr$fields, cols)
    expect_identical(hdr$report_created, "Thu May 28 13:30:51 EDT 2026")
})

test_that("the header is found whatever its length", {
    # A hard-coded line count is what this replaces: the bundled sample
    # and a real download have preambles of different lengths.
    cols <- .ctd_cols()
    row <- "Acetaminophen,D000082,,TNF,7124,,Homo sapiens,9606,x,y,1"
    for (n in c(1L, 12L, 40L)) {
        tmp <- .write_fake_ctd(cols, row, n_preamble = n)
        hdr <- ctdR:::.parse_ctd_header(tmp)
        expect_identical(hdr$fields, cols)
        unlink(tmp)
    }
})

test_that("a file with no Report created line parses with NA provenance", {
    tmp <- .write_fake_ctd(.ctd_cols(),
        "Acetaminophen,D000082,,TNF,7124,,Homo sapiens,9606,x,y,1",
        report_created = NA)
    on.exit(unlink(tmp))

    hdr <- ctdR:::.parse_ctd_header(tmp)
    expect_true(is.na(hdr$report_created))
    expect_identical(hdr$fields, .ctd_cols())
})

test_that("import_CTD errors on a file with a header but no records", {
    tmp_file <- .write_fake_ctd(.ctd_cols())
    on.exit(unlink(tmp_file))

    expect_error(import_CTD(tmp_file), "empty or has too few rows")
})

test_that("import_CTD caches data correctly", {
    skip_on_cran()
    skip_if_not_installed("readr")
    skip_if_not_installed("AnnotationDbi")
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("plyr")

    # Create a minimal fake CTD CSV file
    data_rows <- c(
        "Acetaminophen,D000082,,TNF,7124,,Homo sapiens,9606,increases expression,increases expression,12345",
        "Acetaminophen,D000082,,IL6,3569,,Homo sapiens,9606,increases expression,increases expression,12346",
        "Benzene,D001554,,TP53,7157,,Homo sapiens,9606,affects binding,affects binding,12347",
        "Aspirin,D001241,,PTGS2,5743,,Mus musculus,10090,decreases activity,decreases activity,12348"
    )

    tmp_file <- .write_fake_ctd(.ctd_cols(), data_rows)

    # Temporarily redirect the cache to a temp location
    tmp_cache <- file.path(tempdir(), "ctdR_test_cache")
    options(ctdR.cache = tmp_cache)
    on.exit({
        options(ctdR.cache = NULL)
        unlink(tmp_file)
        unlink(tmp_cache, recursive = TRUE)
    })

    expect_message(import_CTD(tmp_file), "Reading CTD")

    bfc <- ctdR:::.ctd_bfc(tmp_cache)
    expect_true(ctdR:::.ctd_cache_has(bfc, "chemicals"))
    expect_true(ctdR:::.ctd_cache_has(bfc, "ChemicalName_GeneEntrezIds"))
    expect_true(ctdR:::.ctd_cache_has(bfc, "ChemicalName_GeneSymbols"))

    # Verify cached data - only human chemicals
    chemicals <- ctdR:::.ctd_cache_load(bfc, "chemicals")
    expect_true("D000082" %in% chemicals$ChemicalID)
    expect_true("D001554" %in% chemicals$ChemicalID)
    expect_false("D001241" %in% chemicals$ChemicalID)

    ChemicalName_GeneEntrezIds <- ctdR:::.ctd_cache_load(
        bfc, "ChemicalName_GeneEntrezIds")
    expect_true("D000082" %in% names(ChemicalName_GeneEntrezIds))
    expect_equal(length(ChemicalName_GeneEntrezIds[["D000082"]]), 2)
})
