#' @title Generate volcano plot
#'
#' @param data
#' @export volcano_plot
#'
#' @import ggplot2 plotly

volcano_plot <- function(input, data, pathway_dic){

  data$circToLin <- NULL
  data$sum_lin <- NULL
  data$sum_junction <- NULL
  cat(colnames(data))
  col_nr <- grep(colnames(data),pattern = gsub(input$volcano_col,pattern = "-",replacement = "_"))[1]
  col_vals <- c(data[,col_nr])
  #req(length(unique(col_vals)) < 100)
  if(sum(sapply(col_vals,is.numeric))>10){
    if(input$log_scale){
      if(min(col_vals,na.rm = T) >= 0&max(col_vals,na.rm = T) <= 1){
        col_vals <- (-1)*log10(col_vals)
        legend_name <- paste0("-log10(",colnames(data)[col_nr],")")
      }
      else{
        col_vals <- sign(col_vals)*log10(abs(col_vals))#if continues variable, then create a broad variation
        legend_name <- paste0("sign*log10(",colnames(data)[col_nr],")")
      }
    }
    else {
      col_vals <- col_vals
      legend_name <- colnames(data)[col_nr]
    }
  }
  else {
    if(length(unique(col_vals)) > 30){ #avoid matching with overwhelming annotation
      showNotification("the matching term exceeds 30 unique identifiers and does not work as color annotation",type = "message")
      col_nr <- grep(colnames(data),pattern = "gene_biotype")[1]
    }
    col_vals <- factor(col_vals)
    legend_name <- colnames(data)[col_nr]
  }
  experiment_id <- gsub(input$experiment_id,pattern = "-",replacement = "_")
  g <- ggplot(data = data) +
    geom_hline(yintercept = -log10(eval(parse(text = input$p))), linetype="dashed") +
    geom_vline(xintercept = -log2(input$fc), linetype="dashed") +
    geom_vline(xintercept = log2(input$fc), linetype="dashed") +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    theme(
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25))) +
    labs(color = legend_name) +
    if(input$experiment_id!=""){
      temp <- colnames(data)[grep(colnames(data),pattern = experiment_id)]
      new_padj <- colnames(data)[grep(colnames(data),pattern = experiment_id)[grep(temp,pattern = "padj")]]
      new_fc <- colnames(data)[grep(colnames(data),pattern = experiment_id)[grep(temp,pattern = "log2")]]
      temp2 <- paste0("`),color = col_vals, text=paste(gene_symbol, col_vals), key = paste(ensembl_gene_id,gene_symbol,gene_biotype,wiki_link,",paste(names(pathway_dic), collapse = ','),",sep = ';')))")
      eval(parse(text = paste0("geom_point(data=data,aes(x=`",new_fc,"`, y=-log10(`",new_padj,temp2)))
    }
    else{
      temp2 <- paste0("),color = col_vals, text=paste(gene_symbol, col_vals), key = paste(ensembl_gene_id,gene_symbol,gene_biotype,wiki_link,",paste(names(pathway_dic), collapse = ','),",sep = ';')))")
      eval(parse(text = paste0("geom_point(data=data,aes(x=log2FoldChange, y=-log10(padj",temp2)))
    }
  if(sum(sapply(col_vals,is.numeric))>10){
    if(!isColor(input$col_high)|!isColor(input$col_low)){
      showNotification("One or both colors are not considered real",type = "message")
      h <- "red"
      l <- "blue"
    }
    else {
      h <- input$col_high
      l <- input$col_low
    }
    #req(isColor(input$col_high)&isColor(input$col_low))
    g <- g + scale_color_gradient(high = h, low = l)
  }
  p <- ggplotly(g,tooltip = "text") %>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "volcanoplot",
        width = 1500,
        height = 1000
      )
    ) %>% event_register('plotly_selected')
  return(p)
}
