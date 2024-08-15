gsea <- function(ChemicalName_GeneSymbols,entrez_ids )
{

  # Create a named numeric vector for fgsea
  stats <- -log10(entrez_ids[,2] + 1e-10)
  names(stats) <- entrez_ids$entrez_ids
  stats <- sort(stats, decreasing = TRUE)

  # Run fgsea
  fgsea_results <- fgsea::fgsea(pathways = ChemicalName_GeneSymbols, stats = stats, nperm = 10000)

  fgsea_results <- as.data.frame(fgsea_results)

  # sort by padj asc
  fgsea_results <- fgsea_results[order(fgsea_results$padj),]

  fgsea_results$foldEnrichment <- abs(fgsea_results$ES) / mean(fgsea_results$ES)

  # change pathway column name to substance
  names(fgsea_results)[1] <- "ChemicalID"

  # Add ChemicalName
  fgsea_results <- merge(fgsea_results, chemicals, by="ChemicalID")

  entrez_ids$Symbol <- rownames(entrez_ids)
  # Add Enriched GENE
  fgsea_results$Enriched_GENE <- sapply(fgsea_results$ChemicalID, function(x) paste(entrez_ids[entrez_ids$entrez_ids %in% ChemicalName_GeneSymbols[[x]],"Symbol"], collapse=", "))

  return(fgsea_results)

}
