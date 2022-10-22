#' Title
#'
#' @param sc_data single-cell RNA-seq expression matrix of related disease.
#'
#' @return imputed single-cell matrix
#' @export
#' @importFrom Rmagic magic

impute_sc.data <- function(sc_data){
  sc_data <- t(sc_data)
  keep_cols <- colSums(sc_data > 0) > 10
  sc_data <- sc_data[,keep_cols]
  imputed_sc_data <- as.matrix(magic(sc_data, genes="all_genes"))
  imputed_sc_data <- t(imputed_sc_data)
  return(imputed_sc_data)
}
