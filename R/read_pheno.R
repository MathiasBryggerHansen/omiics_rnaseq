#' @title Read the file including phenotypes
#'
#' @param pheno the path to the pheontype file
#'
#' @export read_pheno
#'
#' @return a dataframe with the phenotypes

read_pheno <- function(pheno){
  ph <- readLines(pheno)
  temp <- data.frame(matrix(unlist(sapply(ph, FUN =  strsplit,split = "\t")),ncol = 2, byrow = T))
  temp[,1] <- make.names(gsub(temp[,1],pattern = " ",replacement = "."))
  return(temp[-1,])
}
