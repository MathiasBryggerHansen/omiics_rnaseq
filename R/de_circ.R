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


de_circ <- function(input, data, data_lin, pheno, ensembl2id, i){
  sampleNumber = length(pheno[[input$sample_col1]])
  #pheno <- data.frame(pheno[[1]])#data.frame(c(rep("S34F",3),rep("wt",3))
  p <- pheno[[input$sample_col1]]
  ##if the SA/SD method is used, sum up these:
  if(sum(grepl(colnames(data),pattern = "_SA$")) == sampleNumber){
    junctions <- grepl(colnames(data),pattern = "_SD$|_SA$")&!grepl(colnames(data),pattern = "total")
    junction_data <- data[,junctions][,seq(1,length(p),2)] + data[,junctions][,seq(2,length(p),2)]
    colnames(junction_data) <- p
    #junction_data$gene_symbol <- data$gene_symbol
    data$sum_lin <- data$total_SD + data$total_SA
    data$sum_junction <- data$total_junction
    data$gene_symbol <- gsub(data$circRNA_name,pattern = ".*[0-9]_",replacement = "") #this should remove the circRNA tag
    data <- merge(data, ensembl2id, by = "gene_symbol")
    #junction_data <- merge(junction_data, ensembl2id, by = "gene_symbol")
    #junction_data$wikigene_id <- NULL
    #junction_data$gene_biotype <- NULL

  }
  else {#if CIRI2 with BSJ/LIN
    data$ensembl_gene_id <- gsub(data$gene_id,pattern = "\\..*",replacement = "")
    data <- merge(data, ensembl2id, by = "ensembl_gene_id")
    row.names(data) <- paste(data$Internal_circRNA_ID, data$gene_symbol, sep = "_")
    junctions <- grepl(colnames(data),pattern = "CIRI2.circRNAs.txt_BSJ$")
    linear <- grepl(colnames(data),pattern = "CIRI2.circRNAs.txt_LIN$")
    colnames(data) = gsub(pattern = ".CIRI2.circRNAs.txt_BSJ|.CIRI2.circRNAs.txt_LIN", "", colnames(data))
    junction_data <- data[,junctions]
    #junction_data$ensembl_gene_id <- data$ensembl_gene_id
    linear_data <- data[,linear]
    data$sum_lin <- apply(linear_data,1, FUN = sum)
    data$sum_junction <- data$Total_BSJ
  }
  data$circToLin <- data$sum_junction/data$sum_lin
  print(head(junction_data))
  print(head(data_lin))
  row.names(data_lin) <- make.names(gsub("\\..+$", "",row.names(data_lin)),unique = T)
  data <- rbind(junction_data, data_lin)
  res <- count2deseq_analysis(input = input, countdata = junction_data,pheno = pheno, i = i)
  print(head(res[["test"]]))
  res[["test"]] <- res[["test"]][row.names(res)%in%row.names(junction_data),]
  res[["test"]] <- res[["test"]][order(row.names(res[["test"]])),]
  data <- data[order(row.names(data)),] #make sure that the ids match in order from info file and DE res
  res[["circ_info"]] <- data[,c("ensembl_gene_id","circToLin","sum_lin","sum_junction")]
  return(res)
}
