# import_CTD() hands this function one row per interaction, so a chemical
# with several interactions arrives several times under the same name.
# That is the normal shape of CTD, not a data-quality problem, and must
# not raise the warning reserved for an ID carrying two different names.

test_that("repeated rows of the same chemical do not warn", {
    chem <- data.frame(
        ChemicalID = c("D001554", "D001554", "D001554", "D000082"),
        ChemicalName = c("Benzene", "Benzene", "Benzene", "Acetaminophen"),
        stringsAsFactors = FALSE)
    expect_no_warning(out <- .deduplicate_chemicals(chem))
    expect_identical(out$ChemicalID, c("D001554", "D000082"))
})

test_that("an ID with two different names still warns, keeping the first", {
    chem <- data.frame(
        ChemicalID = c("D001554", "D001554", "D001554"),
        ChemicalName = c("Benzene", "Benzene", "Benzol"),
        stringsAsFactors = FALSE)
    expect_warning(out <- .deduplicate_chemicals(chem),
        "1 ChemicalID\\(s\\) appear with more than one ChemicalName")
    expect_identical(out$ChemicalName, "Benzene")
})
