#' Title
#'
#' @param graph PPI data
#' @param distance_matrix matrix calculated by igraph
#' @param cell_signature.list all cell signatures
#' @param phenotye_sample_signature.list signatures of Samples belonging to the same phenotype
#'
#' @return a list of proximity between all cells and the phenotype
#' @export
#'

cells_distance_score <- function(graph, distance_matrix, cell_signature.list, phenotye_sample_signature.list){
  cells_scores <- NULL
  for (i in 1:length(cell_signature.list)) {
    cells_scores[i] <- cell_to_phenotype(graph, distance_matrix, cell_signature.list[[i]], phenotye_sample_signature.list)
  }

  names(cells_scores) <- names(cell_signature.list)
  return(cells_scores)
}
