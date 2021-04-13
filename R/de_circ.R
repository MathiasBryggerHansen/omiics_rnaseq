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
  if(sum(grepl(colnames(data),pattern = "_SA$")&!grepl(colnames(data),pattern = "total")) == sampleNumber){
    print("SA")
    row.names(data) <- data$circRNA_name
    junctions <- grepl(colnames(data),pattern = "_SD$|_SA$")&!grepl(colnames(data),pattern = "total")
    junction_data <- data[,junctions][,seq(1,length(p)*2,2)] + data[,junctions][,seq(2,length(p)*2,2)]
    colnames(junction_data) <- p
    data$sum_lin <- data$total_SD + data$total_SA
    data$sum_junction <- data$total_junction
    data$gene_symbol <- gsub(data$circRNA_name,pattern = ".*[0-9]_",replacement = "") #this should remove the circRNA tag
    data <- merge(data, ensembl2id, by = "gene_symbol", all.x =T)
    print("gene_symbol")
    print(head(data))
    #circ2ensembl <- data[,c("ensembl_gene_id","circRNA_name")]
  }
  else {#if CIRI2 with BSJ/LIN
    data$ensembl_gene_id <- gsub(data$gene_id,pattern = "\\..*",replacement = "")
    data <- merge(data, ensembl2id, by = "ensembl_gene_id", all.x =T)
    data$circRNA_name <- paste(data$Internal_circRNA_ID, data$gene_symbol, sep = "_")
    row.names(data) <- data$circRNA_name #paste(data$Internal_circRNA_ID, data$gene_symbol, sep = "_")
    junctions <- grepl(colnames(data),pattern = "CIRI2.circRNAs.txt_BSJ$")
    linear <- grepl(colnames(data),pattern = "CIRI2.circRNAs.txt_LIN$")
    colnames(data) = gsub(pattern = ".CIRI2.circRNAs.txt_BSJ|.CIRI2.circRNAs.txt_LIN", "", colnames(data))
    junction_data <- data[,junctions]
    linear_data <- data[,linear]
    data$sum_lin <- apply(linear_data,1, FUN = sum)
    data$sum_junction <- data$Total_BSJ
  }
  circ2ensembl <- data[,c("ensembl_gene_id","circRNA_name")]
  data$circToLin <- data$sum_junction/data$sum_lin
  row.names(data_lin) <- make.names(gsub("\\..+$", "",row.names(data_lin)),unique = T)
  circIDs <- row.names(junction_data)
  linIDs <- row.names(data_lin)
  all_data <- rbind(junction_data, data_lin)
  print(head(all_data))
  row.names(all_data) <- c(circIDs, linIDs)
  res <- count2deseq_analysis(input = input, countdata = all_data,pheno = pheno, i = i)
  all_IDs <- row.names(res[["test"]])
  res[["test"]] <- data.frame(res[["test"]])
  res[["test"]] <- res[["test"]][all_IDs%in%circIDs,]
  res[["test"]]$circRNA_name <- row.names(res[["test"]])
  res[["test"]] <- merge(res[["test"]],circ2ensembl,by = "circRNA_name")
  res[["test"]] <- res[["test"]][order(row.names(res[["test"]])),]
  data <- data[order(row.names(data)),] #make sure that the ids match in order from info file and DE res
  res[["circ_info"]] <- data[,c("circRNA_name","ensembl_gene_id","circToLin","sum_lin","sum_junction")]
  return(res)
}
