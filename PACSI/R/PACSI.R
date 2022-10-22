


#' Title
#'
#' binary group indicator vector,
#' @param bulk_data Bulk expression matrix of related disease. Each row represents a gene and each column represents a sample.
#' @param phenotye_labels Phenotype annotation of each bulk sample. It should be a binary variable.
#' @param sc_data single-cell RNA-seq expression matrix of related disease. Each row represents a gene and each column represents a sample.
#' @param ppi_data Network data that fits the input of the igraph package
#' @param sc_VariableFeatures The number of highly variable genes that need to be retained during single-cell data preprocessing.
#' @param times Number of perturbations.
#' @param ncores Number of CPU cores for parallel operations.
#' @param sample_signature_threshold The size of the sample signature.
#' @param cell_signature_threshold The size of the cell signature.
#' @param corr_p_threshold The P-value threshold for calculating the adjacency matrix.
#' @param cell_population_FDR_threshold The threshold for predicting cells with phenotype 1.
#' @param spatial If the input of sc_data is a seurat object for the spatial transcriptome, let spatial =T.
#' @param impute If you want to fill in missing values in single-cell data, let impute =T.
#' @param bulk_VariableFeatures  If you need to extract the highly variable genes of the bulk matrix, let bulk_VariableFeatures be the number of features such as bulk_VariableFeatures = 2000.
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom igraph distances
#' @importFrom igraph V
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom stats p.adjust
#' @return asdfas
#' @export


PACSI <- function(bulk_data, phenotye_labels, sc_data, ppi_data = NULL, sc_VariableFeatures = 2000, times = 100, ncores = 4,spatial=F,impute=F,bulk_VariableFeatures = F,
                  sample_signature_threshold = 150, cell_signature_threshold = 150,corr_p_threshold=0.05,
                  cell_population_FDR_threshold = 0.05){


  print("build the graph")
  graph <- Extract_max_connected_subgraph(ppi_data, sc_data,corr_p_threshold=corr_p_threshold)

  print("Data preprocessing")
  if(spatial==T){
    sc_data <- as.matrix(sc_data@assays$SCT@scale.data);gc()
  }else{
    if(impute==T){
      imputed_sc_data <- suppressMessages(impute_sc.data(sc_data))
    }
    sc_data <- suppressWarnings(Seurat::CreateSeuratObject(as.matrix(sc_data)))
    sc_data <- Seurat::FindVariableFeatures(object = sc_data, verbose = F, nfeatures=sc_VariableFeatures)
    sc_data <- Seurat::ScaleData(object = sc_data, verbose = F)
    sc_data <- as.matrix(sc_data@assays$RNA@scale.data);gc()
  }
  bulk_data <- bulk_data[apply(bulk_data, 1, function(x) sum(x > 0) > 0.5*ncol(bulk_data)),]
  if(bulk_VariableFeatures != F){
    bulk_variable <- Seurat::FindVariableFeatures(object = bulk_data, verbose = F)
    bulk_variable <- bulk_variable[order(bulk_variable$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
    bulk_variable <- rownames(bulk_variable)[1:bulk_VariableFeatures]
    bulk_data <- bulk_data[match(bulk_variable,rownames(bulk_data)),]
  }
  bulk_data <- scale(bulk_data)

  print("Calculate the cell signatures and sample signatures")
  sample_ranked.matrix <- calculate_ranked.matrix(bulk_data)
  sample_signature.matrix <- calculate_signature(sample_ranked.matrix, sample_signature_threshold)
  sample_signature.list <- list()
  sample_signature.list <- apply(sample_signature.matrix, 2, function(i){
    tmp <- rownames(sample_signature.matrix)[which(i == 1)]
    return(tmp)
  })
  phenotye_sample_signature.list <- extract_phenotye_sample_signature(sample_signature.list, phenotye_labels)

  cell_ranked.matrix <- calculate_ranked.matrix(sc_data)
  cell_signature.matrix <- calculate_signature(cell_ranked.matrix, cell_signature_threshold)
  cell_signature.list <- list()
  cell_signature.list <- apply(cell_signature.matrix, 2, function(i){
    tmp <- rownames(cell_signature.matrix)[which(i == 1)]
    return(tmp)
  })

  if(length(ppi_data) == 0){
    cell_uniprot_signature.list_screened <- signature.list_screened(cell_signature.list, graph)
    phenotye_sample_uniprot_signature.list_screened <- signature.list_screened(phenotye_sample_signature.list, graph)
  }else{
    cell_uniprot_signature.list  <- suppressMessages(symbol_to_uniprot(cell_signature.list))
    cell_uniprot_signature.list_screened <- signature.list_screened(cell_uniprot_signature.list, graph)
    phenotye_sample_uniprot_signature.list  <- suppressMessages(symbol_to_uniprot(phenotye_sample_signature.list))
    phenotye_sample_uniprot_signature.list_screened <- signature.list_screened(phenotye_sample_uniprot_signature.list, graph)
  }
  rm("bulk_data","cell_ranked.matrix","cell_signature.list","cell_signature.matrix","phenotye_sample_signature.list","sample_ranked.matrix","sample_signature.list","sample_signature.matrix");gc()



  print("Calculate the distance between cells and the phenotype")
  distance_matrix <- as.matrix(igraph::distances(graph, V(graph), V(graph)))
  cells_scores <- NULL;
  cl <- makeCluster(ncores);
  clusterEvalQ(cl, library(igraph));
  clusterExport(cl,c("cell_to_sample","cell_to_phenotype"));
  cells_scores <- parSapply(cl, cell_uniprot_signature.list_screened, function(i, graph, distance_matrix, phenotye_sample_uniprot_signature.list_screened){
    tmp <- cell_to_phenotype(graph, distance_matrix, i, phenotye_sample_uniprot_signature.list_screened)
    return(tmp)
  }, graph, distance_matrix, phenotye_sample_uniprot_signature.list_screened);
  stopCluster(cl);
  names(cells_scores) <- names(cell_uniprot_signature.list_screened)



  print("permutation")
  permuted_cells_scores_matrix <- NULL
  cl <- makeCluster(ncores);
  clusterEvalQ(cl, library(igraph));
  clusterExport(cl,c("permute_cell_signature", "cells_distance_score", "cell_to_sample","cell_to_phenotype"));
  system.time({
    permuted_cells_scores_matrix <- parSapply(cl,1:times, function(i, graph, distance_matrix, cell_uniprot_signature.list_screened, phenotye_sample_uniprot_signature.list_screened){
      tmp_permuted_cell_signature.list <- permute_cell_signature(graph, cell_uniprot_signature.list_screened)
      tmp_permuted_cells_scores <- cells_distance_score(graph, distance_matrix, tmp_permuted_cell_signature.list, phenotye_sample_uniprot_signature.list_screened)
      return(tmp_permuted_cells_scores)
    },graph, distance_matrix, cell_uniprot_signature.list_screened, phenotye_sample_uniprot_signature.list_screened)
  })
  stopCluster(cl)


  print("calculate p value")
  p.value <- calculate_pvalue(cells_scores, permuted_cells_scores_matrix, times)
  FDR <- p.adjust(p.value, method = "BH")
  cell_phenotype_labels <- cal_cell_phenotype_labels(FDR, cell_population_FDR_threshold)



  parameter <- list( ppi_data = ppi_data, permutation.times = times, ncores = ncores,sc_VariableFeatures= sc_VariableFeatures,spatial=spatial,impute=impute,bulk_VariableFeatures = bulk_VariableFeatures,
                     sample_signature_threshold = sample_signature_threshold, cell_signature_threshold = cell_signature_threshold,
                     adjacent.matrix_p_threshold=corr_p_threshold, cell_population_FDR_threshold = cell_population_FDR_threshold)
  result <- data.frame(cells_scores = cells_scores, p.value = p.value, FDR=FDR, cell_phenotype_labels = cell_phenotype_labels)
  return(list(parameter=parameter, sc_data=sc_data, graph=graph,result=result,cell_uniprot_signature.list_screened=cell_uniprot_signature.list_screened,phenotye_sample_uniprot_signature.list_screened=phenotye_sample_uniprot_signature.list_screened,permuted_cells_scores_matrix=permuted_cells_scores_matrix))
}
