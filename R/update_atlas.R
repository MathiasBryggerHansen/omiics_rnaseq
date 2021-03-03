#' @title Append data from the Atlas db to dataset
#'
#' @param gene_data the dataset
#' @param atlas_id
#' @param input
#'
#' @export update_atlas

update_atlas <- function(atlas_id = F, gene_data, input){
  atlas <- search_Atlas(atlas_id, input = input)
  if(length(names(atlas))>200){
    showNotification(paste("Found",length(names(atlas)),"results, try a more specific term"),type = "message")
  }
  req(length(names(atlas))<200)
  if(length(atlas)==0){
    return(gene_data)
  }
  fir <- T
  for(n in names(atlas)){
    temp <- atlas[[n]]
    temp <- temp[,c("padj","log2FoldChange")] #change, since there may be more

    colnames(temp) <- paste(colnames(temp),n,sep = "_")
    temp$ensembl_gene_id <- row.names(temp)
    if(fir){
      temp2 <- merge(gene_data, temp, by = "ensembl_gene_id",all.x = T)
      fir <- F
    }
    else{
      temp2 <- merge(temp2, temp, by = "ensembl_gene_id",all.x = T)
    }
  }
  return(temp2)
}
