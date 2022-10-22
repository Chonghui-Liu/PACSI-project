

#' Title
#'
#' @param uniprot_signature.list uniprot signature
#' @param graph network
#' @importFrom igraph V
#' @return IDs that exists on the network.
#' @export
#'

signature.list_screened <- function(uniprot_signature.list, graph){
  uniprot_signature.list_screened <- list()
  for (i in 1:length(uniprot_signature.list)) {
    uniprot_signature.list_screened[[i]] <- uniprot_signature.list[[i]][is.element(uniprot_signature.list[[i]],as.character(names(V(graph))))]
  }
  names(uniprot_signature.list_screened) <- names(uniprot_signature.list)
  return(uniprot_signature.list_screened)
}
