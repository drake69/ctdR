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

test_that("import_CTD errors on wrong file format (missing columns)", {
    # Create a CSV with wrong columns (e.g. a different CTD file)
    header_lines <- paste0("# line ", seq_len(27))
    col_names <- c("WrongCol1", "WrongCol2", "WrongCol3")
    dummy_row <- paste(rep("dummy", length(col_names)), collapse = ",")
    data_rows <- "val1,val2,val3"

    tmp_file <- tempfile(fileext = ".csv")
    writeLines(c(header_lines, paste(col_names, collapse = ","), dummy_row, data_rows), tmp_file)
    on.exit(unlink(tmp_file))

    expect_error(
        import_CTD(tmp_file),
        "Not a valid CTD file"
    )
    expect_error(
        import_CTD(tmp_file),
        "Missing columns"
    )
})

test_that("import_CTD errors on empty file", {
    header_lines <- paste0("# line ", seq_len(27))
    col_names <- c(
        "ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID", "Interaction",
        "InteractionActions", "PubMedIDs"
    )

    tmp_file <- tempfile(fileext = ".csv")
    writeLines(c(header_lines, paste(col_names, collapse = ",")), tmp_file)
    on.exit(unlink(tmp_file))

    expect_error(
        import_CTD(tmp_file),
        "empty or has too few rows"
    )
})

test_that("import_CTD errors on wrong file format (missing columns)", {
  # Create a CSV with wrong columns (e.g. a different CTD file)
  header_lines <- paste0("# line ", seq_len(27))
  col_names <- c("WrongCol1", "WrongCol2", "WrongCol3")
  dummy_row <- paste(rep("dummy", length(col_names)), collapse = ",")
  data_rows <- "val1,val2,val3"

  tmp_file <- tempfile(fileext = ".csv")
  writeLines(c(header_lines, paste(col_names, collapse = ","), dummy_row, data_rows), tmp_file)
  on.exit(unlink(tmp_file))

  expect_error(
    import_CTD(tmp_file),
    "does not match the expected CTD"
  )
  expect_error(
    import_CTD(tmp_file),
    "Missing columns"
  )
})

test_that("import_CTD errors on empty file", {
  header_lines <- paste0("# line ", seq_len(27))
  col_names <- c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
                 "GeneForms", "Organism", "OrganismID", "Interaction",
                 "InteractionActions", "PubMedIDs")

  tmp_file <- tempfile(fileext = ".csv")
  writeLines(c(header_lines, paste(col_names, collapse = ",")), tmp_file)
  on.exit(unlink(tmp_file))

  expect_error(
    import_CTD(tmp_file),
    "empty or does not contain enough data"
  )
})

test_that("import_CTD caches data correctly", {
    skip_on_cran()
    skip_if_not_installed("readr")
    skip_if_not_installed("AnnotationDbi")
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("plyr")

    # Create a minimal fake CTD CSV file
    header_lines <- paste0("# line ", seq_len(27))
    col_names <- c(
        "ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID", "Interaction",
        "InteractionActions", "PubMedIDs"
    )
    dummy_row <- paste(rep("dummy", length(col_names)), collapse = ",")
    data_rows <- c(
        "Acetaminophen,D000082,,TNF,7124,,Homo sapiens,9606,increases expression,increases expression,12345",
        "Acetaminophen,D000082,,IL6,3569,,Homo sapiens,9606,increases expression,increases expression,12346",
        "Benzene,D001554,,TP53,7157,,Homo sapiens,9606,affects binding,affects binding,12347",
        "Aspirin,D001241,,PTGS2,5743,,Mus musculus,10090,decreases activity,decreases activity,12348"
    )

    tmp_file <- tempfile(fileext = ".csv")
    writeLines(c(header_lines, paste(col_names, collapse = ","), dummy_row, data_rows), tmp_file)

    # Temporarily override user_cache_dir to use a temp location
    tmp_cache <- file.path(tempdir(), "ctdR_test_cache")
    original_cache_fn <- rappdirs::user_cache_dir
    assignInNamespace("user_cache_dir", function(...) tmp_cache, ns = "rappdirs")
    on.exit({
        assignInNamespace("user_cache_dir", original_cache_fn, ns = "rappdirs")
        unlink(tmp_file)
        unlink(tmp_cache, recursive = TRUE)
    })

    expect_message(import_CTD(tmp_file), "Reading CTD")

    expect_true(file.exists(file.path(tmp_cache, "chemicals.rda")))
    expect_true(file.exists(file.path(tmp_cache, "ChemicalName_GeneEntrezIds.rda")))
    expect_true(file.exists(file.path(tmp_cache, "ChemicalName_GeneSymbols.rda")))

    # Verify cached data - only human chemicals
    load(file.path(tmp_cache, "chemicals.rda"))
    expect_true("D000082" %in% chemicals$ChemicalID)
    expect_true("D001554" %in% chemicals$ChemicalID)
    expect_false("D001241" %in% chemicals$ChemicalID)

    load(file.path(tmp_cache, "ChemicalName_GeneEntrezIds.rda"))
    expect_true("D000082" %in% names(ChemicalName_GeneEntrezIds))
    expect_equal(length(ChemicalName_GeneEntrezIds[["D000082"]]), 2)
})
