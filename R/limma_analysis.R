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

limma_analysis <- function(countdata_norm,phenotypes){#expects normalized data #,p = "control|normal|reference" #,clin = NULL
  res <- list()
  rownames(countdata_norm) <- make.names(gsub("\\..+$", "",row.names(countdata_norm)),unique = T)#gsub("\\..+$", "",rownames(countdata_norm))
  phenotypes <- gsub(phenotypes,pattern = "-",replacement = "_")
  if(length(unique(phenotypes))==2&sum(grepl(x = phenotypes,pattern = "control|normal|reference",ignore.case = T))>0){
    ref <- names(table(phenotypes))[grepl(names(table(phenotypes)),pattern = "control|normal|reference",ignore.case = T)]
    p_id <- names(table(phenotypes))[!grepl(names(table(phenotypes)),pattern = p,ignore.case = T)] #phenotype
    design <- cbind(Control=1,CasevsControl=ifelse(phenotypes==ref,1,0) )
    fit <- lmFit(countdata, design)
    fit <- eBayes(fit, trend=TRUE)
    temp <- topTable(fit, coef=ncol(design),number = nrow(countdata))
  }
  else{
    f <- factor(phenotypes, levels=unique(phenotypes))
    design <- model.matrix(~0+f)
    fit <- lmFit(countdata_norm, design)
    a <- combinations(phenotypes, n = length(unique(phenotypes)),r = 2)
    a <- paste(paste0("f",a[,1]),paste0("f",a[,2]),sep = "-")
    c <- paste(a,collapse = ", ")
    d <- "contrast.matrix <- makeContrasts("
    e <- ",levels=design)"
    eval(parse(text=paste(d,c,e,collapse = "")))
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    temp <- topTable(fit2,number = nrow(countdata_norm))#, coef=ncol(design),number = nrow(countdata))
  }
  if("ID"%in%colnames(temp)){#failsafe, this occurs with duplicate ids, but this should be corrected
    row.names(temp) <- make.names(temp$ID,unique = T)
  }
  if(length(a) == 1){
    temp$log2FoldChange <- temp$logFC
    res[["test"]] <- temp[,c("log2FoldChange","AveExpr","P.Value","adj.P.Val","B")]
    colnames(res[["test"]]) <- c("log2FoldChange","baseMean","pvalue","padj","B")
    res[["test"]] <- res[["test"]][,c("padj","log2FoldChange","baseMean","pvalue","B")]
  }
  else{
    temp$log2FoldChange <- apply(temp[,1:length(a)],1,FUN = function(r) {r[which.max(abs(r))]})
    logFC_ids <- paste0("log2FoldChange_",colnames(temp[,1:length(a)]))
    res[["test"]] <- temp[,c("log2FoldChange",colnames(temp[,1:length(a)]),"AveExpr","P.Value","adj.P.Val","F")]
    colnames(res[["test"]]) <- c("log2FoldChange",logFC_ids ,"baseMean","pvalue","padj","F")
    res[["test"]] <- res[["test"]][,c("padj","log2FoldChange",logFC_ids,"baseMean","pvalue","F")]
  }
  countdata_norm <- countdata_norm[order(row.names(countdata_norm)),]
  res[["norm_counts"]] <- countdata_norm
  res[["phenotypes"]] <- phenotypes
  return(res)
}
