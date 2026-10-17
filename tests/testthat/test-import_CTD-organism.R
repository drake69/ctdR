# import_CTD(organism = ) keeps the species asked for and nothing else.
# The default must reproduce the human-only import of every earlier
# version, and the choice must reach the provenance record, otherwise two
# caches built from the same file look identical.

.org_fixture <- function() {
    cols <- c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID", "Interaction",
        "InteractionActions", "PubMedIDs")
    rows <- c(
        "Benzene,D001554,,TP53,7157,,Homo sapiens,9606,x,affects^binding,1",
        "Benzene,D001554,,CYP1A1,1543,,Homo sapiens,9606,x,increases^expression,2",
        "Benzene,D001554,,Cyp1a1,13076,,Mus musculus,10090,x,increases^expression,3",
        "Benzene,D001554,,Cyp1a1,24296,,Rattus norvegicus,10116,x,increases^expression,4",
        "Benzene,D001554,,Ahr,11622,,,,x,affects^binding,5"
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

.with_org_cache <- function(code) {
    tmp_cache <- file.path(tempdir(), paste0("ctdR_org_", basename(tempfile())))
    options(ctdR.cache = tmp_cache)
    on.exit({
        .restore_ctd_cache()
        unlink(tmp_cache, recursive = TRUE)
    }, add = TRUE)
    force(code)
}

.cached_entrez <- function() {
    sets <- .ctd_cache_load(.ctd_bfc(), "ChemicalName_GeneEntrezIds")
    sort(unlist(sets, use.names = FALSE))
}

test_that(".check_organism accepts numbers, strings and NULL", {
    expect_null(.check_organism(NULL))
    expect_identical(.check_organism(9606), "9606")
    expect_identical(.check_organism(c(9606, 10090, 9606)), c("9606", "10090"))
    expect_identical(.check_organism("10116"), "10116")
})

test_that(".check_organism rejects values that are not taxonomy IDs", {
    expect_error(.check_organism(character(0)), "NCBI Taxonomy IDs")
    expect_error(.check_organism(NA), "NCBI Taxonomy IDs")
    expect_error(.check_organism("Homo sapiens"), "not valid: Homo sapiens")
    expect_error(.check_organism(9606.5), "not valid")
    expect_error(.check_organism(list(9606)), "NCBI Taxonomy IDs")
})

test_that("the default keeps only human interactions, as before", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        expect_message(import_CTD(f), "Filtered to 2 human interactions")
        expect_identical(.cached_entrez(), c("1543", "7157"))
        expect_identical(ctd_provenance()$organism, "9606")
    })
})

test_that("a vector of species keeps exactly those species", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        expect_message(import_CTD(f, organism = c(9606, 10090)),
            "No orthology mapping")
        expect_identical(.cached_entrez(), c("13076", "1543", "7157"))
        expect_identical(ctd_provenance()$organism, c("9606", "10090"))
    })
})

test_that("NULL keeps every species, including unrecorded ones", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        expect_message(import_CTD(f, organism = NULL), "species not recorded 1")
        expect_identical(.cached_entrez(),
            c("11622", "13076", "1543", "24296", "7157"))
        prov <- ctd_provenance()
        expect_true(is.na(prov$organism))
        expect_output(print(prov), "Organism: +all species")
    })
})

test_that("symbol sets keep only genes with a human Entrez ID", {
    # No orthology: the mouse and rat Cyp1a1 must not turn into CYP1A1.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        suppressMessages(import_CTD(f, organism = NULL))
        syms <- .ctd_cache_load(.ctd_bfc(), "ChemicalName_GeneSymbols")
        expect_setequal(syms$gene, c("TP53", "CYP1A1"))
    })
})

test_that("a species absent from the file warns, and none at all stops", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        expect_warning(suppressMessages(
            import_CTD(f, organism = c(9606, 7955))),
            "No interactions in the file for organism 7955")
        expect_error(suppressMessages(import_CTD(f, organism = 7955)),
            "No interactions in the file for organism 7955")
    })
})

test_that("a record written before the field existed prints as human", {
    expect_identical(.format_organism(NULL),
        "9606 (Homo sapiens; not recorded)")
    expect_identical(.format_organism(c("9606", "10090")), "9606, 10090")
})
