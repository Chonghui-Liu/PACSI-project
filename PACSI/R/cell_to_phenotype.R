#' Title
#'
#' @param graph PPI network
#' @param distance_matrix matrix calculated by igraph
#' @param cell_signature cell signature
#' @param phenotye_sample_signature.list signatures of Samples belonging to the same phenotype.
#' @importFrom stats na.omit
#' @return proximity between a cell and the phenotype.
#' @export
#'

cell_to_phenotype <- function(graph, distance_matrix,cell_signature,phenotye_sample_signature.list){
  cell_to_sample_vector <- NULL
  for (i in 1:length(phenotye_sample_signature.list)) {
    cell_to_sample_vector[i] <- cell_to_sample(graph, distance_matrix,cell_signature,phenotye_sample_signature.list[[i]])
  }
  cell_to_phenotype_distance <- mean(na.omit(cell_to_sample_vector))
  return(cell_to_phenotype_distance)
}
