#' Filter DE results on log2foldchange and p-value, and includes 10 percent of "non-significant"
#'
#' @param data the dateframe
#' @export filter_results
#' @return the filtered dataframe

filter_results <- function(input, data){
  data_sign <- data[data$padj < eval(parse(text = input$p)) & abs(data$log2FoldChange) > log2(input$fc),]
  print(str(data_sign))
  data_not_sign <- data[data$padj > eval(parse(text = input$p)) | abs(data$log2FoldChange) < log2(input$fc),]
  print(head(data_not_sign))
  data_not_sign_10 <- data_not_sign[sample(seq(1,nrow(data_not_sign)),size = round(nrow(data_not_sign)/10)),]
  return(rbind(data_sign,data_not_sign_10))
}
