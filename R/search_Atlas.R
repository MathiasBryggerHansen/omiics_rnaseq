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

search_Atlas <- function(input, atlas_id = F, probe_library = probe_library()) {#c("padj","log2FoldChange","baseMean","pvalue","t","B") from limma_analysis
  res <- list()
  p <- NULL
  #uses probe library to convert probes to ensembl gene id
  try({
    if(atlas_id){
      acc <- strsplit(gsub(input$Atlas_ids,pattern = " ",replacement = ""),split = ",")
      allExps <- getAtlasData(acc)
    }
    else{
      showNotification(paste("Contacting server... This may take a while"),type = "message")
      atlasRes <- searchAtlasExperiments(properties = input$atlas_query, species = input$species)
      allExps <- getAtlasData(atlasRes$Accession)
      showNotification(paste("Found",length(names(allExps)),"Atlas search results..."),type = "message")
    }

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
  })

  return(res)
}
