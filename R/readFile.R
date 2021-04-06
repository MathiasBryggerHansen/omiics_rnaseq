readFile <- function(f) {
  pheno_data <- tryCatch({
  read.table(f,header = T,stringsAsFactors = F)},#read_pheno(input[[paste0("phenotype",d)]][["datapath"]]) #phenotype data is always assumed to be tabulated, the function handles some errors in read.csv
  warning = function(cond) {
    showNotification("The count file could not be read returning:",type = "message")
    showNotification(cond,type = "message")
    return(NULL)},
  error = function(cond) {
    showNotification("The count file could not be read returning:",type = "message")
    showNotification(cond,type = "message")
    return(NULL)}
  )
  return(pheno_data)
}
