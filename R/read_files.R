#' function used by input_d
#' if files are not combined this reads from a folder containing the single samples
#'
#' @param dir the directory of the folder
#' @param d the delimiter
#' @export read_files
#'
#' @return the combined dataset

read_files <- function(input,dir){
  for (f in list.files(dir)){
    #f <- paste(dir,f, sep = "/")
    if(grepl(f,pattern = ".gz",fixed = T)){
      #untar(zipfile = f) #paste0(dir,"/",
      temp <- read.csv(gzfile(paste0(dir,"/",f)), sep = input[[paste0("sep",d)]], header = T,stringsAsFactors=FALSE)
      #f <- strsplit(f,split = ".gz")[[1]]
    }
    else if(grepl(f,pattern = ".zip",fixed = T)){
      #unzip(zipfile = f)
      # f <- strsplit(f,split = ".zip")[[1]]
      temp <- read.csv(unz(paste0(dir,"/",f)), sep = input[[paste0("sep",d)]], header = T,stringsAsFactors=FALSE)
    }
    n <- ncol(temp)
    temp <- temp[grepl(temp[,1],pattern = "ENS|MST"),] #remove non-genes
    temp <- temp[order(temp[,1]),]
    temp[,1] <- gsub("\\..+$", "",temp[,1])
    temp <- temp[!duplicated(temp[,1]),] #removes PAR_Y

    if(exists(x = "res")){
      res <- cbind(res,as.numeric(temp[[n]]))
    }
    else {
      res <- as.numeric(temp[,n])
    }
  }
  res <- data.frame(res)
  row.names(res) <- temp[,1]
  samples <- gsub(pattern = "\\..+$",replacement = "",list.files(dir))
  colnames(res) <- samples
  return(res)
}
