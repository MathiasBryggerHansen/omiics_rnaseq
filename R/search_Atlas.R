#' @title Search the Atlas database
#'
#' @param input from UI, inclludes species, atlas id, atlas query
#' @param probe_library a library translating different ids
#'
#' @import limma ExpressionAtlas edgeR
#'
#' @export search_Atlas
#'
#' @return data.frame of DE analysis

search_Atlas <- function(input, atlas_id = F, probe_library = probe_library()) {
  res <- list()
  p <- NULL
  #uses probe library to convert probes to ensembl gene id
  if(atlas_id){
    acc <- unlist(strsplit(gsub(input$Atlas_ids,pattern = " ",replacement = ""),split = ","))
    print(acc)
    allExps <- tryCatch({
      getAtlasData(acc)
    },
    error=function(cond){
      showNotification("No results found for the given query",type = "message")
      return(NULL)
    })
  }
  else{
    showNotification(paste("Contacting server... This may take a while"),type = "message")
    allExps <- tryCatch({
      getAtlasData(acc)
    },
    error=function(cond){
      showNotification("No results found for the given query",type = "message")
      return(NULL)
    })
  }
  if(is.null(allExps)){
    return(NULL)
  }
  for(fi in names(allExps)){
    try({
      showNotification(paste("Analyzing",fi),type = "message")
      type <- names(allExps[[fi]])
      temp <- allExps[[fi]][[1]]
      pheno <- temp$AtlasAssayGroup #this gets very messy -> gsub(paste0(temp$AtlasAssayGroup,temp$clinical_information,temp$disease),pattern = " ",replacement = ".") #AtlasAssayGroup should be a general ID
      if(length(unique(pheno))==1|length(unique(pheno))>8|length(pheno)==length(unique(pheno))){
        next
      }
      if(grepl(type,"rnaseq")){#Data from atlas is either raw RNAseq counts or normalized array data
        d0 <- DGEList(assays(temp)$counts)
        if(input$gene_filter){
          d0 <- d0[filterByExpr(d0, group=pheno),, keep.lib.sizes=FALSE]
        }
        d0 <- calcNormFactors(d0, method = "TMM")# #TMM normalization #assays(counts) - only from atlas
        f <- factor(pheno, levels=unique(pheno))

        mm <- model.matrix(~0+f)
        temp <- voom(d0, mm)
      }
      else{
        temp <- exprs(temp)
        print("her")
      }
      p_id <- gsub(fi,pattern = "-",replacement = "_")
      ids <- row.names(temp)
      if(!grepl(ids[1],pattern = "ENS")){
        ids <- probe_library$ensembl_gene_id[match(x = ids, probe_library$probe)]
      }
      row.names(temp) <- make.names(ids,unique = T)
      print("w3")
      if(sum(is.na(ids)) == length(ids)){ #One reason for this error is that not all the R files (from ebi) are structured correctly with wrong feature names. Some of this could be fixed by using the .txt tables instead.
        showNotification("Gene IDs of experiment not recognized", type = "message")
        stop("Gene IDs of experiment not recognized")
      }
      res[[p_id]] <- limma_analysis(countdata = temp,phenotypes = pheno, control = input$control1, auto = T)$test
    })
  }
  ##clean up
  remove(allExps)
  return(res)
}
