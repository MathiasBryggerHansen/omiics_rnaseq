#' Installs and loads all packages needed

install_load <- function(packages){
  for (p in packages) {
    if (p %in% rownames(installed.packages())) {
      library(p, character.only=TRUE)
    } else {
      BiocManager::install(p)
      library(p,character.only = TRUE)
    }
  }
}
