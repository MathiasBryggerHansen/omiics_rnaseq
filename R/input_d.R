#' Reads in the count files
#' Handles if the files are separate for each sample, the zipped folder is unzipped and the data combined.
#'
#' @param input UI data including the amount datasets (nfiles)
#' @export input_d
#'
#' @return list of "count", "pheno" and "circ" (if circRNA data is included)

input_d <- function(input){
  files <- list()
  for( d in 1:input$nfiles){
    print("her")

    if(input[[paste0("combined",d)]]){#test for combined text file
      h <- readLines(input[[paste0("count",d)]][["datapath"]], n = 1)
      print("her2")
      if(grepl(h, pattern = "featureCounts")){ #if the dataset is raw from featureCounts
        print("her3")
        counts_data <- read.csv2(input[[paste0("count",d)]][["datapath"]], sep = input[[paste0("sep",d)]], header = T,comment.char = "#",skip = 1,stringsAsFactors = F)
        print(head(counts_data))
        colnames(counts_data) <- sapply(colnames(counts_data),FUN = function(x) strsplit(x, split = ".", fixed = T)[[1]][1]) #set sample ids
        counts_data <- counts_data[,-c(seq(2,6))] #remove Chr	Start	End	Strand	Length
      }
      else {
        counts_data <- read.csv2(input[[paste0("count",d)]][["datapath"]], sep = input[[paste0("sep",d)]], header = T,comment.char = "!",stringsAsFactors = F) #comment.char = "!" in CEL files
      }
    }
    else{
      d_path <- input[[paste0("count",d)]][["datapath"]]
      temp_dir <- paste0("./temp",d)
      unlink(temp_dir, recursive=TRUE,force = T) #remove if exists
      dir.create(temp_dir)
      unzip(d_path,exdir = temp_dir,overwrite = T)
      dir_name <- list.files(temp_dir)
      counts_data <- read_files(paste0(temp_dir,"/",dir_name),d)
    }
    if(!is.null(input[[paste0("circRNA",d)]])){#test circ data
      circ_data <- read.csv2(input[[paste0("circRNA",d)]][["datapath"]], sep = input[[paste0("sep",d)]], header = T,comment.char = "!",stringsAsFactors = F) #comment.char = "!" in CEL files
    }
    else{
      circ_data <- NULL
    }

    if(input$gene_filter){
      counts_data <- counts_data[apply(X = counts_data,1, function(x) var(x)!=0),] #remove zero variance genes, Warning in var(x) : NAs introduced by coercion
    }
    print("4")
    pheno_data <- read.table(input[[paste0("phenotype",d)]][["datapath"]],header = T,stringsAsFactors = F)#read_pheno(input[[paste0("phenotype",d)]][["datapath"]]) #phenotype data is always assumed to be tabulated, the function handles some errors in read.csv
    print(pheno_data)
    if(input$gene_id_col & input[[paste0("combined",d)]]){
      row.names(counts_data) <- make.names(counts_data[,1],unique = T)
      counts_data[,1] <- NULL
    }
    pheno_keep <- grepl(colnames(counts_data),pattern = paste(pheno_data[[input$sample_col]],collapse = "|"))
    ids <- row.names(counts_data)
    if(!grepl(ids[1],pattern = "ENS")){
      ids <- probe_library()$ensembl_gene_id[match(x = ids, probe_library()$probe)]
      row.names(counts_data) <- make.names(ids,unique = T)
    }
    counts_data <- counts_data[,c(pheno_keep)]#colnames(counts_data)%in%pheno_data[[2]]] #remove samples not in pheno
    colnames(counts_data) <- pheno_data[[input$sample_col]]
    counts_data <- counts_data[order(row.names(counts_data)),] #a way to secure compatible order????
    files[[paste0("count",d)]] <- counts_data
    files[[paste0("pheno",d)]] <- pheno_data
    files[[paste0("circRNA",d)]] <- circ_data
  }
  return(files)
}
