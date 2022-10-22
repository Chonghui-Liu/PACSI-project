#' Title
#'
#' @param sc_data single-cell RNA-seq expression matrix of related disease. Each row represents a gene and each column represents a sample.
#' @param corr_p_threshold The P-value threshold for calculating the adjacency matrix.
#'
#' @return Adjacency matrix for single-cell data
#' @importFrom Hmisc rcorr
#' @export
#'
#' @examples test=cal_adj_matrix(case1_sc_data_1to100,0.05)
cal_adj_matrix <- function(sc_data, corr_p_threshold){
  adj_matrix <- rcorr(t(as.matrix(sc_data)))
  adj_matrix <- adj_matrix$P
  adj_matrix[adj_matrix < corr_p_threshold] <- 1
  adj_matrix[adj_matrix != 1] <- 0
  adj_matrix[is.na(adj_matrix)] <- 0
  return(adj_matrix)
}
