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


readUrl <- function(url) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      message("This is the 'try' part")

      readLines(con=url, warn=FALSE)
      # The return value of `readLines()` is the actual value
      # that will be returned in case there is no condition
      # (e.g. warning or error).
      # You don't need to state the return value via `return()` as code
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", url))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      message(paste("Processed URL:", url))
      message("Some other message at the end")
    }
  )
  return(out)
}

search_Atlas <- function(input, atlas_id = F, probe_library = probe_library()) {#c("padj","log2FoldChange","baseMean","pvalue","t","B") from limma_analysis
  res <- list()
  p <- NULL
  #uses probe library to convert probes to ensembl gene id
    if(atlas_id){
      acc <- strsplit(gsub(input$Atlas_ids,pattern = " ",replacement = ""),split = ",")
      allExps <- tryCatch({
        print("trying")
        getAtlasData(acc)
      },
      error=function(cond){
        print("error")
        showNotification("No results found for the given query",type = "message")
        return(NULL)
      })
    }
    else{
      showNotification(paste("Contacting server... This may take a while"),type = "message")
      # atlasRes <- searchAtlasExperiments(properties = input$atlas_query, species = input$species)
      # tryCatch({
      #   print("trying")
      #   allExps <- getAtlasData(acc)
      # },
      # error=function(cond){
      #   print("error")
      #   showNotification("No results found for the given query",type = "message")
      #   return(NULL)
      # })
      allExps <- tryCatch({
        print("trying")
        getAtlasData(acc)
      },
      error=function(cond){
        print("error")
        showNotification("No results found for the given query",type = "message")
        return(NULL)
      })
    }
  if(is.null(allExps)){
    return(NULL)
  }
    print("continuing")
    for(fi in names(allExps)){
      try({
        showNotification(paste("Analyzing",fi),type = "message")
        #p <- "control|normal|reference"
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
          #counts <- cpm(counts,log = T)
          f <- factor(pheno, levels=unique(pheno))#phenotypes[[1]]

          mm <- model.matrix(~0+f)
          temp <- voom(d0, mm)
        }
        else{
          temp <- exprs(temp)
        }
        p_id <- gsub(fi,pattern = "-",replacement = "_")
        ids <- row.names(temp)
        if(!grepl(ids[1],pattern = "ENS")){
          ids <- probe_library$ensembl_gene_id[match(x = ids, probe_library$probe)]
        }
        row.names(temp) <- make.names(ids,unique = T)
        res[[p_id]] <- limma_analysis(countdata = temp,phenotypes = pheno)$test
      })
    }
    ##clean up
    remove(allExps)
  return(res)
}
