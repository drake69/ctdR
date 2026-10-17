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

test_that(".check_organism accepts numbers, strings and \"all\"", {
    expect_identical(.check_organism("all"), "all")
    expect_identical(.check_organism("ALL"), "all")
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
    expect_error(.check_organism(NULL), "\"all\"")
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
            "not an orthology mapping")
        expect_identical(.cached_entrez(), c("13076", "1543", "7157"))
        expect_identical(ctd_provenance()$organism, c("9606", "10090"))
    })
})

test_that("\"all\" keeps every species, including unrecorded ones", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        expect_message(import_CTD(f, organism = "all"), "species not recorded 1")
        expect_identical(.cached_entrez(),
            c("11622", "13076", "1543", "24296", "7157"))
        prov <- ctd_provenance()
        expect_identical(prov$organism, "all")
        expect_output(print(prov), "Organism: +all species")
    })
})

test_that("symbol sets pool species by gene symbol", {
    # Mouse and rat Cyp1a1 count as CYP1A1, mouse Ahr as AHR: a match by
    # name, as on the CTD website. Entrez sets keep each species' gene.
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        suppressMessages(import_CTD(f, organism = "all"))
        syms <- .ctd_cache_load(.ctd_bfc(), "ChemicalName_GeneSymbols")
        expect_setequal(syms$gene, c("TP53", "CYP1A1", "AHR"))
        expect_identical(sum(syms$gene == "CYP1A1"), 1L)
        expect_identical(.cached_entrez(),
            c("11622", "13076", "1543", "24296", "7157"))
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

test_that("a human-only import keeps the symbols it always had", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        suppressMessages(import_CTD(f))
        syms <- .ctd_cache_load(.ctd_bfc(), "ChemicalName_GeneSymbols")
        expect_setequal(syms$gene, c("TP53", "CYP1A1"))
    })
})

test_that("the interaction-type filter pools symbols the same way", {
    skip_on_cran()
    skip_if_not_installed("org.Hs.eg.db")
    f <- .org_fixture()
    on.exit(unlink(f))

    .with_org_cache({
        suppressMessages(import_CTD(f, organism = "all"))
        gs <- suppressMessages(
            .filter_gene_sets(.ctd_cache_dir(), "increases^expression"))
        expect_setequal(gs$symbols$gene, "CYP1A1")
    })
})

test_that(".gene_symbol_map keeps the human-only rule unless pooling", {
    skip_if_not_installed("org.Hs.eg.db")
    d <- data.frame(GeneID = c("7157", "999999999", "18671"),
        GeneSymbol = c("TP53", "OLDNAME", "Abcb1a"),
        OrganismID = c("9606", "9606", "9606"))
    # Human-only: org.Hs.eg.db alone, unknown IDs dropped, as before.
    m <- .gene_symbol_map(d)
    expect_identical(m[["7157"]], "TP53")
    expect_true(is.na(m[["999999999"]]))
    expect_true(is.na(m[["18671"]]))
    # Pooling: decided on the Entrez ID, not on the row's organism, so a
    # mouse gene recorded on a human interaction keeps its CTD symbol.
    p <- .gene_symbol_map(d, pool_species = TRUE)
    expect_identical(p[["7157"]], "TP53")
    expect_identical(p[["18671"]], "ABCB1A")
})

test_that("the species message lists at most ten species", {
    many <- data.frame(Organism = paste("sp", 1:12), OrganismID = as.character(1:12))
    expect_message(.filter_organism(many, "all"), "and 2 more")
})
