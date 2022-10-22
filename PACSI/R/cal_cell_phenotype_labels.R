#' Title
#'
#' @param p.value The p-value of the perturbation test
#' @param p.value_threshold The threshold for predicting cells with phenotype 1.
#'
#' @return Predicted cell phenotypic labels
#' @export
#'

cal_cell_phenotype_labels <- function(p.value, p.value_threshold){
  cell_phenotype_labels <- NULL
  for (i in 1:length(p.value)) {
    if(p.value[i] <= p.value_threshold){
      cell_phenotype_labels[i] = 1
    }else{
      cell_phenotype_labels[i] = 0
    }
  }
  return(cell_phenotype_labels)
}
