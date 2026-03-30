ora <- function(ChemicalName_GeneSymbols,entrez_ids){

  gene_list <- entrez_ids
  gene2pathway <- ChemicalName_GeneSymbols

  # Utilizza l'enricher per l'analisi ORA
  ora_results <- clusterProfiler::enricher(gene = gene_list, TERM2GENE = gene2pathway)

  # browser()
  ora_results <- ora_results@result

  # rename p.adjust as padj
  colnames(ora_results)[colnames(ora_results)=="p.adjust"] <- "padj"

  # calculate foldEnrichment using BgRatio and GeneRatio
  # dplyr::mutate(ora_results, foldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
  ora_results$foldEnrichment <- DOSE::parse_ratio(ora_results$GeneRatio) / DOSE::parse_ratio(ora_results$BgRatio)

  return(ora_results)
}
