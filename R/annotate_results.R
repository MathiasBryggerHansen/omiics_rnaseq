#' @title Annotation of DE results
#'
#' @description Translates Ensembl gene IDs to gene symbol, information on Transcription factor status and the pathways in which the gene is involved. It adds secondary DE analysis results along with annotation from chosen files.
#'
#' @param data a DE result data.frame including either ensembl_gene_id as a column (for circRNA data) or row names for linRNA.
#'
#' @param input UI input, includes species
#'
#' @param circ boolean indicating if the data is of circRNA data structure or linear.
#'
#' @param ensembl2id a reactive data.frame including gene symbol information
#'
#' @param pathway_dic a reactive data.frame with the pathways involved for each gene
#'
#' @export annotate_results
#'
#' @return data.frame

annotate_results <- function(input, data,ensembl2id,pathway_dic,circ = F){
  if(!circ){#is it circRNA data?
    data$ensembl_gene_id <- row.names(data)
  }
  print("annotate_merge")
  data <- merge(data,ensembl2id,by = "ensembl_gene_id",all.x = T)#ensembl2id() includes : gene_symbol, wiki_id, biotype
  print("merge_done")
  data$gene_symbol <- gsub(make.names(ifelse(is.na(data$gene_symbol)|data$gene_symbol=="",data$ensembl_gene_id,data$gene_symbol),unique = T),pattern = ".",replacement = "-",fixed = T)
  if(input$species == "human"|input$species == "mouse"){ #Transcription factor annotation
    tf_data <- read.table(paste0("data/trrust_rawdata.",input$species,".tsv"),header = T) #this should be more general
    tf_data$position <- NULL
    data$TF <- data$gene_symbol%in%tf_data$gene_symbol
    print("TF")
  }
  data <- append_secondary_data(input, data,de_results = gene_results_de())#append the secondary datasets
  data <- append_annotation(input, data)#append annotation
  if(input$species == "human"|input$species == "mouse"){ #add pathway information
    paths <- data.frame(matrix(rep(NA,length(data$gene_symbol)*length(names(pathway_dic))),ncol = length(names(pathway_dic))))
    colnames(paths) <- names(pathway_dic)
    print("dict")
    for(db in 1:length(names(pathway_dic))){
      for(s in 1:length(data$gene_symbol)){
        paths[s,db] <- paste(pathway_dic[[names(pathway_dic)[db]]][[data$gene_symbol[s]]],collapse = ", ")
      }
    }
    for(db in names(pathway_dic)){
      data[[db]] <- paths[,db]
    }
  }
  else{
    for(db in names(pathway_dic)){
      data[[db]] <- "-"
    }
  }
  wiki_link <- sapply(data$wikigene_id, function(x) ifelse(is.na(x),"-",paste0(c("https://www.wikigenes.org/e/gene/e/",x,".html"),collapse = "")))
  refs <- paste0("<a href='",  wiki_link, "' target='_blank'>",wiki_link,"</a>")
  data$wiki_link <- refs
  data$wikigene_id <- NULL
  data <- data[order(data$padj),]
  print("annotated")
  print(head(data))
  return(data)
}
