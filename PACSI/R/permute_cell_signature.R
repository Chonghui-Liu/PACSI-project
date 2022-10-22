#' Title
#'
#' @param graph PPI data
#' @param cell_uniprot_signature.list all cell signatures
#' @importFrom igraph get.vertex.attribute
#' @return random cell signatures
#' @export
#'
#'
permute_cell_signature <-  function(graph, cell_uniprot_signature.list){
  permuted_cell_signature.list <- list()
  for (i in 1:length(cell_uniprot_signature.list)) {
    permuted_cell_signature.list[[i]] <- sample(get.vertex.attribute(graph)[[1]], size=length(cell_uniprot_signature.list[[i]]), replace = TRUE)
  }
  names(permuted_cell_signature.list) <- names(cell_uniprot_signature.list)
  return(permuted_cell_signature.list)
}
