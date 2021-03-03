#' @title Differential Expression analysis on circRNA
#'
#' @description Structures the data, and calls the count2deseq_analysis function
#'
#' @param data count data
#' @param pheno phenotype dataframe
#'
#' @return list of results
#'
#' @export de_circ


de_circ <- function(data, pheno){
  sampleNumber = 6#length(pheno[[1]])
  pheno <- data.frame(c(rep("S34F",3),rep("wt",3)))
  data$ensembl_gene_id <- gsub(data$gene_id,pattern = "\\..*",replacement = "")
  data$gene_id <- NULL
  data <- merge(data,ensembl2id(),all.x = T, by = "ensembl_gene_id")
  colnames(data) = gsub(pattern = ".CIRI2.circRNAs.txt_BSJ|.CIRI2.circRNAs.txt_LIN", "", colnames(data))
  row.names(data) <- paste(data$gene_symbol,data$Internal_circRNA_ID,sep = "_")
  counts <- data[,c(14:(2*sampleNumber + 13))] #(data[,c(13:(2*sampleNumber + 12))])
  j_index = seq(1,(2*sampleNumber),2)
  lin_index = seq(2,(2*sampleNumber),2)
  rawdata <- counts[,j_index]
  counts$ensembl_gene_id <- data$ensembl_gene_id
  counts$sum_lin <- rowSums(counts[lin_index])
  counts$sum_junction <- rowSums(counts[j_index])#######overwrite pheno
  counts$circToLin <- counts$sum_junction/counts$sum_lin #2*counts[j_index]/(2*counts[j_index]+counts[lin_index])
  res <- count2deseq_analysis(countdata = rawdata,pheno = pheno)
  res[["circ_info"]] <- counts[,c("ensembl_gene_id","circToLin","sum_lin","sum_junction")]
  return(res)
}
