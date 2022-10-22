#' Title
#'
#' @param graph PPI data
#' @param distance_matrix matrix calculated by igraph
#' @param cell_signature cell signature
#' @param sample_signature sample signature
#' @importFrom igraph V
#' @return proximity between a cell and a sample
#' @export
#'
#'
cell_to_sample <- function(graph, distance_matrix, cell_signature, sample_signature){
  if(length(cell_signature) == 0 | length(sample_signature) == 0){
    cell_to_sample_distance <- NA
  }else{
    cell_to_sample_distance <- mean(apply(as.matrix(distance_matrix[which(is.element(names(V(graph)), cell_signature)),which(is.element(names(V(graph)), sample_signature))]), 1, min))
  }
  return(cell_to_sample_distance)
}
