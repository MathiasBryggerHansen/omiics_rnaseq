#' @title Appending secondary DE results to primary data set.
#'
#' @param data primary data set on which to append secondary DE results
#' @param de_results a reactive data frame including the secondary DE results
#' @param ensembl_lift a text file to translate ensembl ids between species
#'
#' @export append_secondary_data
#'
#' @return data.frame

append_secondary_data <- function(input, data,de_results = gene_results_de(),ensembl_lift = "human_mouse_ensembl.txt"){#This is used both when selecting from Volcanoplot and when generating gene_data$df, all.x = T is used otherwise it would annotate all genes
  if(input$nfiles == 1){
    return(data)
  }
  else{
    if(input$nfiles > 1){
      for(i in 2:input$nfiles){
        temp <- de_results[[toString(i)]][["test"]][,c("padj","log2FoldChange")]
        colnames(temp) <- c(paste0("padj_",input[[paste0("phen",i)]]),paste0("log2FoldChange_",input[[paste0("phen",i)]]))
        temp$ensembl_gene_id <- row.names(temp)
        if(sum(row.names(data)%in%temp$ensembl_gene_id)<10){
          humanToMouse <- read.table(ensembl_lift,header = T)
          if(grepl(data$ensembl_gene_id[1],pattern = "ENSG")){
            temp <- merge(temp,humanToMouse, by.x = "ensembl_gene_id",by.y = "ensembl_gene_id_mouse")
            temp$ensembl_gene_id <- temp$ensembl_gene_id_human
            temp$ensembl_gene_id_human <- NULL
            temp$ensembl_gene_id_mouse <- NULL
          }
          else{
            temp <- merge(temp,humanToMouse, by.x = "ensembl_gene_id",by.y = "ensembl_gene_id_human")
            temp$ensembl_gene_id <- temp$ensembl_gene_id_mouse
            temp$ensembl_gene_id_human <- NULL
            temp$ensembl_gene_id_mouse <- NULL
          }
        }
        data <- merge(data,temp,by = "ensembl_gene_id",all.x = T)
      }
    }
  }
  return(data)
}
