#' Connects to the Atlas API, and tests for overlapping expression
#'
#' @param A vector of ensembl gene ids to test against the database
#'
#' @return data.frame including dataset id and link
#'
#' @export testGenes

testGenes <- function(gene_vector, ebi_id, input) {
  out <- tryCatch(
    {
      gene_vector <- as.vector(na.omit(gsub("\\..+$", "",x = gene_vector))) #remove ENS.X subscripts
      out <- read.csv(url(paste0('https://www.ebi.ac.uk/fg/gsa/api/tsv/getOverlappingComparisons/',ebi_id[[input$species]],"/",paste(gene_vector,collapse = "%20"))),sep = "\t",skip = 2)
      out <- out[,c("EXPERIMENT","ADJUSTED.P.VALUE","EXPERIMENT_URL")]
      colnames(out) <- c("dataset","padj","link")
      out <- out[out$padj < 0.05,]
      out <- out[order(out$padj,decreasing = F),]
      out$link <- paste0("<a href='",  out$link, "' target='_blank'>",out$link,"</a>")
      return(out)

    },
    error=function(cond) {
      message(paste("There may be too few significant genes:", paste0('https://www.ebi.ac.uk/fg/gsa/api/tsv/getOverlappingComparisons/',ebi_id[[input$species]],"/",paste(gene_vector,collapse = "%20"))))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", paste0('https://www.ebi.ac.uk/fg/gsa/api/tsv/getOverlappingComparisons/',ebi_id[[input$species]],"/",paste(gene_vector,collapse = "%20"))))
      message("Here's the original warning message:")
      message(cond)
      return(NULL)
    }
  )
  return(out)
}
