#' Title
#'
#' @param sample_signature.list all sample signatures
#' @param phenotye_labels pheontype labels
#'
#' @return signatures of Samples belonging to the same phenotype
#' @export
#'

extract_phenotye_sample_signature <- function(sample_signature.list, phenotye_labels){
  phenotye_sample_signature.list <- sample_signature.list[which(phenotye_labels == 1)]
  return(phenotye_sample_signature.list)
}
