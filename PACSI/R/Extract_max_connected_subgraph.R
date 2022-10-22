#' Title
#'
#' @param ppi_data Network data
#' @param sc_data single-cell RNA-seq expression matrix
#' @param corr_p_threshold he P-value threshold for calculating the adjacency matrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph groups
#' @importFrom igraph induced_subgraph
#' @importFrom igraph is.connected
#' @importFrom igraph components
#' @importFrom dplyr distinct
#' @importFrom utils read.delim
#' @return the largest connected subgraph of the input network data
#' @export
#'

Extract_max_connected_subgraph <- function(ppi_data, sc_data, corr_p_threshold){
  if(length(ppi_data) == 0){
    ppi <- cal_adj_matrix(sc_data,corr_p_threshold=corr_p_threshold)
    g <- graph_from_adjacency_matrix(ppi)
  }else{
    ppi <- read.delim(ppi_data,header=T,sep = ",")
    ppi <- ppi[,c(1,2)]
    gc()

    ppi <- distinct(ppi)
    tem_id <- NULL
    for (i in 1:nrow(ppi)) {
      if(ppi[i,1] == ppi[i,2]){
        tem_id[i] <- 1
      }else{
        tem_id[i] <- 0
      }
    }
    table(tem_id)
    ppi <- ppi[-which(tem_id==1),]

    edge<-data.frame(from= ppi[ ,1],to= ppi[ ,2])
    rm(ppi);gc()
    g<-graph_from_data_frame(edge,directed = FALSE)
  }

  if(is.connected(g)){
    graph <- g
  }else{
    g_components <- groups(components(g))
    graph <- induced_subgraph(g, as.character(g_components[[which(lengths(g_components)==max(lengths(g_components)))]]), impl = c("auto"))
  }

  return(graph)
}
