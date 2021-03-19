#' @title DEseq2 analysis of raw count data
#' @description Includes option of batch correction
#'
#' @param countdata an N x M data.frame with N genes and M samples
#' @param pheno an A x B phenotype data.frame with A samples, column 2 will include batch if present.
#'
#' @import DESeq2
#'
#' @return named list, including the variables "test", "norm_counts", "dds", "phenotypes"
#'
#' @export count2deseq_analysis

count2deseq_analysis <- function(input, countdata, pheno){
  countdata <- as.matrix(countdata)
  res <- list()
  row.names(countdata) <- make.names(gsub("\\..+$", "",row.names(countdata)),unique = T)
  line = gsub("_.$", "",colnames(countdata))
  line = factor(line)
  control <- input[[paste0("control",1)]] #this needs to be adjusted if there are multiple files?
  phenotypes <- factor(pheno[[input$group_col1]])
  print(phenotypes)
  if(!control%in%phenotypes){
    showNotification(paste0("Your control group must match one group ID (",paste(phenotypes, collapse = ", "),")"),type = "message")
  }
  req(control%in%phenotypes)
  #case <- input$case
  if(input$batch_correction&ncol(pheno)>1){#values need to be updated if batch correction is chosen
    batch <- pheno[[input$batch_col1]]
    samples <- data.frame(row.names=colnames(countdata),
                          line=line,
                          phenotypes=phenotypes,
                          batch=batch)

    dds <- DESeq2::DESeqDataSetFromMatrix(countData=countdata, samples, design=~batch + phenotypes)
  }
  else {
    samples <- data.frame(row.names=colnames(countdata),
                          line=line,
                          phenotypes=phenotypes)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=countdata, samples, design=~phenotypes)
  }
  # if(sum(grepl(x = phenotypes,pattern = "control|normal|reference|wt",ignore.case = T))>0){
  #   control <- unique(phenotypes)[grepl(x = unique(phenotypes),pattern = "control|normal|reference|wt",ignore.case = T)]
  #   dds$phenotypes <- relevel(dds$phenotypes, control) #sets the control group
  # }
  dds$phenotypes <- relevel(dds$phenotypes, control) #sets the control group
  dds <- DESeq2::DESeq(dds)
  cases <- unique(phenotypes)
  cases <- factor(cases[cases!=control])
  for(case in unique(cases)){
    test <- DESeq2::results(dds,contrast = c("phenotypes",case,control))
    if(!exists("results")){
      if(length(case) == 1){
        results <- test
        colnames(results) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
        results <- results[,c("baseMean","log2FoldChange","padj")]
        results$ensembl_gene_id <- row.names(results)
      }
      else{
        results <- test[,c("baseMean","log2FoldChange","padj")]
        colnames(results) <- c("baseMean",paste0("log2FoldChange_",case),paste0("padj_",case))
        results$ensembl_gene_id <- row.names(results)
      }
    }
    else {
      colnames(test) <- c("baseMean",paste0("log2FoldChange_",case),"lfcSE","stat","pvalue",paste0("padj_",case))
      test <- test[,c(paste0("log2FoldChange_",case),paste0("padj_",case))]
      test$ensembl_gene_id <- row.names(test)
      results <- merge(results, test, by = "ensembl_gene_id")
    }
  }

  all_res$ensembl_gene_id <- NULL
  res[["test"]] <- results
  res[["norm_counts"]] <- assay(varianceStabilizingTransformation(dds))
  res[["dds"]] <- dds
  res[["phenotypes"]] <- phenotypes
  #all possible combinations of phenotype interactions
  # if(length(unique(phenotypes))>2){
  #   for(p in unique(phenotypes)){##needs to be tested with raw multivariate data!!
  #     for(p2 in unique(phenotypes)){
  #       if(p!=p2&!paste0(p,p2)%in%comparisons_rev){
  #         comparisons_rev <- c(comparisons_rev,paste0(p2,p))
  #         comparisons <- paste0(p,p2)
  #         if(!exists("all_res")){
  #           all_res <- data.frame(results(dds,contrast = c("phenotypes",p,p2)))
  #           colnames(all_res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  #           all_res$ensembl_gene_id <- row.names(all_res)
  #         }
  #         else{
  #           temp <- data.frame(results(dds,contrast = c("phenotypes",p,p2)))
  #           colnames(temp) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  #           colnames(temp) <- paste(colnames(temp),comparisons,sep = "_")
  #           temp$ensembl_gene_id <- row.names(temp)
  #           merge(all_res,temp,by = "ensembl_gene_id")
  #         }
  #       }
  #     }
  #   }
  # }
  # else{
  #   all_res <- data.frame(results(dds))
  #   colnames(all_res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  # }
  return(res)
}
