#' Title
#'
#' @param signature.list symbol ID
#' @importFrom clusterProfiler bitr
#' @return uniprot id
#' @export

symbol_to_uniprot <- function(signature.list){
  uniprot_signature.list <- list()
  for (i in 1:length(signature.list)) {
    uniprot_signature.list[[i]] = tryCatch({
      bitr(signature.list[[i]], fromType = "SYMBOL", toType = "UNIPROT",OrgDb = "org.Hs.eg.db")$UNIPROT
    }, error = function(e) {
      as.character()
    })
  }
  names(uniprot_signature.list) <- names(signature.list)
  return(uniprot_signature.list)
}
