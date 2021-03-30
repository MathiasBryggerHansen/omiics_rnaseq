#' Sets up a fitting regression model and runs a limma analysis
#'
#' @param countdata_norm normalized count data
#' @param phenotypes a vector of phenotype ids
#'
#'
#' @import limma gtools
#'
#' @return list of results
#'
#' @export limma_analysis

limma_analysis <- function(countdata_norm, phenotypes, auto = F, control){#expects normalized data #,p = "control|normal|reference" #,clin = NULL
  res <- list()
  rownames(countdata_norm) <- make.names(gsub("\\..+$", "",row.names(countdata_norm)),unique = T)#gsub("\\..+$", "",rownames(countdata_norm))
  phenotypes <- gsub(phenotypes,pattern = "-",replacement = "_")
  f <- factor(phenotypes, levels=unique(phenotypes))
  a <- "contrast.matrix <- makeContrasts("
  c <- ",levels=design)"
  if(length(unique(f)) == 2){
    if(auto&sum(grepl(x = phenotypes,pattern = "control|normal|reference",ignore.case = T))==1){#is auto, but there is a control (the control parameter is not used in auto)
      control <- names(table(phenotypes))[grepl(names(table(phenotypes)),pattern = "control|normal|reference",ignore.case = T)]
    }
    else if(auto){
      control <- unique(f)[1]
    }
    design <- cbind(Control=1,CasevsControl=ifelse(phenotypes==control,1,0))
    fit <- lmFit(countdata_norm, design)
    fit <- eBayes(fit, trend=TRUE)
    combined <- topTable(fit, coef=ncol(design),number = nrow(countdata_norm))
    combined <- combined[,c("logFC","adj.P.Val","AveExpr")]
    colnames(combined) <- c("log2FoldChange","padj","baseMean")
  }
  else{
    design <- model.matrix(~0+f)
    fit <- lmFit(countdata_norm, design)
    if(auto&sum(grepl(x = phenotypes,pattern = "control|normal|reference",ignore.case = T))==1){#is auto, but there is a control (the control parameter is not used in auto)
      control <- names(table(phenotypes))[grepl(names(table(phenotypes)),pattern = "control|normal|reference",ignore.case = T)]
    }
    else if(auto&sum(grepl(x = phenotypes,pattern = "control|normal|reference",ignore.case = T))!=1){#if there is no clear control and is auto){
      combi <- combinations(phenotypes, n = length(unique(phenotypes)), r = 2) #then all combinations are tested
      combi <- paste(paste0("f",combi[,1]),paste0("f",combi[,2]),sep = "-")
      b <- paste(combi,collapse = ", ")
      cases <- combi #since all combinations are tested, each combination is seen as a "case"
      eval(parse(text=paste(a,b,c,collapse = "")))
    }
    else {
      cases <- unique(phenotypes[phenotypes!=control])
      b <- paste(paste0("f",control), paste0("f",cases), sep = "-")
      b <- paste(b, collapse = ", ")
      eval(parse(text=paste(a,b,c,collapse = "")))
    }
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)
    first <- T
    for(i in 1:length(cases)){
      if(first){
        combined <- topTable(fit = fit, coef = i,number = nrow(countdata_norm))
        combined <- combined[,c("logFC","adj.P.Val","AveExpr")]
        colnames(combined) <- c("log2FoldChange","padj","baseMean")
        combined$ensembl_gene_id <- row.names(combined)
        first <- F
      }
      else{
        temp <- topTable(fit = fit, coef = i,number = nrow(countdata_norm))
        temp <- temp[,c("logFC","adj.P.Val")]
        colnames(temp) <- c(paste("log2FoldChange",i,sep = "_"),paste("padj",i,sep = "_") )
        temp$ensembl_gene_id <- row.names(temp)
        #print(head(temp))
        #print("combined:")
        #print(head(combined))
        combined <- merge(combined,temp,by = "ensembl_gene_id")
      }
    }
    if(auto){#if auto settings only the most significant data is used for each gene
      log2Combined <- combined[,grep(colnames(combined),pattern = "log2Fold")]
      padjCombined <- combined[,grep(colnames(combined),pattern = "padj")]
      log2Top <- apply(log2Combined,1,FUN = function(r) {which.max(abs(r))})
      #log2 <- apply(log2Combined,1,FUN = function(r) {log2Combined[which.max(abs(r))]})
      padjTop <- padjCombined[cbind(1:nrow(padjCombined),log2Top)]
      log2Top <- log2Combined[cbind(1:nrow(log2Combined),log2Top)]
      combined <- data.frame(cbind(log2Top, padjTop, combined$ensembl_gene_id, combined$baseMean))
      colnames(combined) <- c("log2FoldChange","padj","ensembl_gene_id","baseMean")
      combined$log2Top <- as.numeric(combined$log2FoldChange)
      combined$padjTop <- as.numeric(combined$padj)
      combined$baseMean <- as.numeric(combined$baseMean)
    }
  }

  row.names(combined) <- combined$ensembl_gene_id
  combined$ensembl_gene_id <- NULL
  res[["fit"]] <- fit
  res[["test"]] <- combined
  countdata_norm <- countdata_norm[order(row.names(countdata_norm)),]
  res[["norm_counts"]] <- countdata_norm
  res[["phenotypes"]] <- phenotypes
  return(res)
}
