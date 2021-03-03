#' @export isColor

isColor <- function(x){
  res <- try(col2rgb(x),silent=TRUE)
  return(!"try-error"%in%class(res))
}
