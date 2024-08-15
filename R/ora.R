ora <- function(ChemicalName_GeneSymbols,entrez_ids){

  gene_list <- entrez_ids
  gene2pathway <- ChemicalName_GeneSymbols

  # Utilizza l'enricher per l'analisi ORA
  ora_results <- clusterProfiler::enricher(gene = gene_list, TERM2GENE = gene2pathway)

  ora_results <- ora_results@result

  return(ora_results)
}
