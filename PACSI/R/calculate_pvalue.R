#' Title
#'
#' @param cells_scores the proximity between cells and the phenotype of interest.
#' @param permuted_cells_scores_matrix perturbation results.
#' @param times Number of perturbations.
#'
#' @return p value of proximity.
#' @export
#'

calculate_pvalue <- function(cells_scores, permuted_cells_scores_matrix, times){
  cells_scores_pvalue <- NULL
  for (i in 1:nrow(permuted_cells_scores_matrix)) {
    if(is.nan(cells_scores[i])){
      cells_scores_pvalue[i] <- 1
    }else{
      cells_scores_pvalue[i] <- length(which(permuted_cells_scores_matrix[i,] < cells_scores[i])) / times
    }
  }
  return(cells_scores_pvalue)
}
