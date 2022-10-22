#' Title
#'
#' @param ranked.matrix rank-based matrix
#' @param signature_threshold the size of signature
#'
#' @return gene signature
#' @export
#'

calculate_signature <- function(ranked.matrix,signature_threshold){
  signature.matrix<- data.frame()
  signature.matrix <- apply(ranked.matrix, 2, function(i){
    tmp <- NULL
    tmp <- sort(i,decreasing = T)
    tmp <- as.numeric(i >= tmp[signature_threshold])
    return(tmp)
  })
  rownames(signature.matrix) <- rownames(ranked.matrix)
  colnames(signature.matrix) <- colnames(ranked.matrix)
  return(signature.matrix)
}
