#' @title Add annotation from extra annotation file
#'
#' @param d a DE result data.frame including either ensembl_gene_id as a column.
#'
#' @return data.frame
#'
#' @export append_annotation

append_annotation <- function(input, d){
  if(input$afiles == 0){
    return(d)
  }
  else{
    for(i in 1:input$afiles){
      ann <- read.csv(input[[paste0("afile",i)]][["datapath"]], sep = input[[paste0("sep_a",i)]], header = T,fill = T) #"\t"
      r <- as.numeric(input[[paste0("anno_range",i)]])
      d <- merge(d,ann[,c(1,r)],by = "ensembl_gene_id", all.x = T)
    }
  }
  d <- d[!duplicated(d$ensembl_gene_id),] #maybe create a list of variables matching instead
  return(d)
}
