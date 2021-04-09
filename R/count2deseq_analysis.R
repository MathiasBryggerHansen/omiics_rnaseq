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

count2deseq_analysis <- function(input, countdata, pheno, i){
  row.names(countdata) <- make.names(gsub("\\..+$", "",row.names(countdata)),unique = T)
  line = gsub("_.$", "",colnames(countdata))
  line = factor(line)
  res <- list()
  control <- input[[paste0("control",i)]]
  case <- input[[paste0("case",i)]]
  phenotypes <- factor(pheno[[input[[paste0("group_col",i)]]]]) #1 --> i
  print(phenotypes)
  print(case)
  print(control)
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
  dds$phenotypes <- relevel(dds$phenotypes, control) #sets the control group
  dds <- DESeq2::DESeq(dds)
  cases <- as.vector(unlist(phenotypes[!phenotypes%in%control]))
  for(i in unique(cases)){
    test <- data.frame(DESeq2::results(dds,contrast = c("phenotypes",i,control)))
    if(!exists("de_res")){
      if(length(unique(cases)) == 1){
        de_res <- test
        colnames(de_res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
        de_res <- de_res[,c("baseMean","log2FoldChange","padj")]
        de_res$ensembl_gene_id <- row.names(de_res)
      }
      if (i == case){
        de_res <- test[,c("baseMean","log2FoldChange","padj")]
        colnames(de_res) <- c("baseMean",paste0("log2FoldChange"),paste0("padj"))
        de_res$ensembl_gene_id <- row.names(de_res)
      }
      else{
        de_res <- test[,c("baseMean","log2FoldChange","padj")]
        colnames(de_res) <- c("baseMean",paste0("log2FoldChange_",i),paste0("padj_",i))
        de_res$ensembl_gene_id <- row.names(de_res)
      }
    }
    else {
      if (i == case){
        colnames(test) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
        test <- test[,c("log2FoldChange","padj")]
        test$ensembl_gene_id <- row.names(test)
        de_res <- merge(de_res, test, by = "ensembl_gene_id")
      }
      else {
        colnames(test) <- c("baseMean",paste0("log2FoldChange_",i),"lfcSE","stat","pvalue",paste0("padj_",i))
        test <- test[,c(paste0("log2FoldChange_",i),paste0("padj_",i))]
        test$ensembl_gene_id <- row.names(test)
        de_res <- merge(de_res, test, by = "ensembl_gene_id")
      }
    }
  }
  row.names(de_res) <- de_res$ensembl_gene_id #why is this needed?
  de_res$ensembl_gene_id <- NULL
  res[["test"]] <- de_res
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
