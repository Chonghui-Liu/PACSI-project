#' Title
#'
#' @param matrix The original expression matrix.
#'
#' @return rank-based matrix.
#' @export
#'

calculate_ranked.matrix <- function(matrix){
  ranked.matrix <- data.frame()
  ranked.matrix <- t(apply(matrix, 1, rank))
  ranked.matrix <- ranked.matrix/ncol(ranked.matrix)
  rownames(ranked.matrix) <- rownames(matrix)
  colnames(ranked.matrix) <- colnames(matrix)
  return(ranked.matrix)
}
