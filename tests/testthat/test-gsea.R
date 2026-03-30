test_that("gsea returns expected columns", {
  skip_on_cran()
  skip_if_not_installed("fgsea")

  # Create minimal pathway list and ranked gene data
  pathways <- list(
    CHEM1 = c(1, 2, 3),
    CHEM2 = c(4, 5, 6)
  )
  chemicals <- data.frame(
    ChemicalID = c("CHEM1", "CHEM2"),
    ChemicalName = c("Chemical One", "Chemical Two")
  )

  set.seed(42)
  entrez_ids <- data.frame(
    entrez_ids = as.character(1:100),
    pvalue = runif(100)
  )

  result <- gsea(pathways, entrez_ids, chemicals)

  expect_true(is.data.frame(result))
  expect_true("ChemicalID" %in% colnames(result))
  expect_true("ChemicalName" %in% colnames(result))
  expect_true("foldEnrichment" %in% colnames(result))
  expect_true("Enriched_GENE" %in% colnames(result))
})

test_that("gsea sorts results by padj", {
  skip_on_cran()
  skip_if_not_installed("fgsea")

  pathways <- list(
    CHEM1 = c(1, 2, 3),
    CHEM2 = c(4, 5, 6)
  )
  chemicals <- data.frame(
    ChemicalID = c("CHEM1", "CHEM2"),
    ChemicalName = c("Chemical One", "Chemical Two")
  )

  set.seed(42)
  entrez_ids <- data.frame(
    entrez_ids = as.character(1:100),
    pvalue = runif(100)
  )

  result <- gsea(pathways, entrez_ids, chemicals)

  if (nrow(result) > 1) {
    expect_true(all(diff(result$padj) >= 0))
  }
})
