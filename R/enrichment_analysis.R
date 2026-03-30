#' @title Enrichment analysis using CTD
#' @description This function performs enrichment analysis using the Comparative Toxicogenomics Database (CTD) data.
#' @param entrez_ids A data frame with two columns: entrez_ids and a numeric value representing the gene expression.
#' @return A data frame with the following columns: ChemicalID, NES, pval, padj, ES, nMoreExtreme, size, leadingEdge, ChemicalName, foldEnrichment, Enriched_GENE.
#' @importFrom rappdirs user_cache_dir
#' @export
enrichment_CTD <- function(entrez_ids, method="ORA")
{
  # Example of getting a cache directory
  cache_dir <- rappdirs::user_cache_dir("ctdR")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  if(!file.exists( file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda")))
  {
    CTD_chem_gene_ixns <- readr::read_csv("~/Downloads/CTD_chem_gene_ixns.csv",  skip = 27)
    CTD_chem_gene_ixns <- CTD_chem_gene_ixns[-1,]

    # rename first columns ChemicalName
    names(CTD_chem_gene_ixns)[1] <- "ChemicalName"

    organisms <- unique(CTD_chem_gene_ixns$OrganismID)
    organisms <- organisms[!is.na(organisms)]

    # Filter only human
    CTD_chem_gene_ixns <- subset(CTD_chem_gene_ixns, CTD_chem_gene_ixns$OrganismID == 9606)

    # rename first column to ChemicalName
    colnames(CTD_chem_gene_ixns)[1] <- "ChemicalName"

    # create a list of ChemicalName and each Chemical a vector of GeneSymbols
    chemicals <- unique(CTD_chem_gene_ixns$ChemicalID)

    # create a list of GeneSymbols for each ChemicalName
    ChemicalName_GeneEntrezIds <- list()
    ChemicalName_GeneSymbols <- data.frame()
    # Create a list of GeneSymbols for each ChemicalName
    browser()
    # for (i in 1:10) {
    for (i in seq_along(chemicals)) {
      # i <- 1
      chemical <- chemicals[i]
      gene_entrezids <-unique(CTD_chem_gene_ixns[CTD_chem_gene_ixns$ChemicalID == chemical,"GeneID"])
      gene_entrezids <- as.data.frame(gene_entrezids)[,1]
      ChemicalName_GeneEntrezIds[[chemical]] <- gene_entrezids

      # transform to a character vector
      gene_entrezids <- as.character(gene_entrezids)
      tryCatch({
        # convert entrezid to symbol
        gene_symbols <- AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = gene_entrezids,
          column = "SYMBOL",
          keytype = "ENTREZID",
          multiVals = "first")
        # ChemicalName_GeneSymbols[[chemical]] <- as.factor(as.vector(gene_symbols))
        tt <- expand.grid(term=chemical, gene=gene_symbols)
        ChemicalName_GeneSymbols <- plyr::rbind.fill(tt,ChemicalName_GeneSymbols)
      }, error = function(e) {
        # next
      })
    }

    chemicals <- unique(CTD_chem_gene_ixns[,c("ChemicalID","ChemicalName")])
    save(chemicals, file = file.path(cache_dir, "chemicals.rda"))
    # write in the cache package folder the set
    save(ChemicalName_GeneEntrezIds, file = file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"))
    save(ChemicalName_GeneSymbols, file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"))
  }
  else
  {
    load(file = file.path(cache_dir, "chemicals.rda"))
  }

  if(method=="ORA")
  {
    load(file = file.path(cache_dir, "ChemicalName_GeneSymbols.rda"))
    gene_symbols <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = entrez_ids$entrez_ids,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first")
    gene_symbols <- gene_symbols[,"SYMBOL"]
    fgsea_results <- ora(ChemicalName_GeneSymbols,gene_symbols)
    # Add ChemicalName
    fgsea_results <- merge(fgsea_results, chemicals, by.y="ChemicalID", by.x="ID")
    # rename ID as ChemicalID
    colnames(fgsea_results)[colnames(fgsea_results)=="ID"] <- "ChemicalID"
    # remove Description column
    fgsea_results <- fgsea_results[,!colnames(fgsea_results) %in% c("Description")]
  }
  else
  {
    load(file = file.path(cache_dir, "ChemicalName_GeneEntrezIds.rda"))
    entrez_ids <- as.data.frame(entrez_ids)
    # Remove NA values
    entrez_ids <- entrez_ids[!is.na(entrez_ids$entrez_ids),]
    fgsea_results <- gsea(ChemicalName_GeneEntrezIds,entrez_ids)
    # Add ChemicalName
    fgsea_results <- merge(fgsea_results, chemicals, by="ChemicalID")
  }

  # browser()

  return(fgsea_results)
}
