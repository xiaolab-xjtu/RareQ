#'
#' @title Obtain more stable clustering result by applying FindRare on shuffled datasets.
#'
#' @description ConsensusRare identifies more stable rare cell types.
#'
#' @usage ConsensusRare(
#'   sc_object,
#'   assay = "RNA",
#'   reduction = "pca",
#'   dims = 1:50,
#'   k.param = 20,
#'   k = 6,
#'   Q_cut = 0.6,
#'   ratio = 0.2,
#'   reps = 30,
#'   user_num =NULL
#' )
#'
#' @param sc_object A seurat object containing assays, reduction, meta data etc.
#' @param assay The Assay object of the input Seurat object for inferring rare cell types. Default is 'RNA'. User can adjust according to data modality.
#' @param reduction Reduction to use as input for shuffling and nearest neighbor searching
#' @param dims Dimensions of reduction to use as input
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param k The size of neighborhood in computing Q values for rare cell detection. k should be between 5 and the total nearest neighbors (e.g., 20) obtained in processing step. The default value is 6.
#' @param Q_cut The default Q value for retaining rare clusters. Q should be between 0.5 and 1. The default value is 0.6.
#' @param ratio The threshold for merging clusters. When the proportion of a cluster's connections with its neighboring clusters among its total connections is larger than this threshold, merge the two cell groups. The default value is 0.2.
#' @param reps The maximum iterations for data shuffling. Default value is 30.
#' @param user_num User defined cluster number.
#' @return Cluster assignments of cells including major and rare cell types.
#' @export
#' @examples
#' ## Starts from counts
#' library(RareQ)
#' library(Seurat)
#'
#' ## Data prepocessing
#' sc_object <- CreateSeuratObject(count=counts, project = "sc_object", min.cells = 3)
#' sc_object$percent.mt <- PercentageFeatureSet(sc_object, pattern = "^MT-")
#' sc_object <- subset(sc_object, percent.mt<20)
#' sc_object <- NormalizeData(sc_object)
#' sc_object <- FindVariableFeatures(sc_object, nfeatures=2000)
#' sc_object <- ScaleData(sc_object, do.center = T)
#' sc_object <- RunPCA(sc_object, npcs=50)
#' sc_object <- FindNeighbors(object = sc_object,
#'                            k.param = 20,
#'                            compute.SNN = F,
#'                            prune.SNN = 0,
#'                            reduction = "pca",
#'                            dims = 1:50,
#'                            force.recalc = F, return.neighbor = T)
#' ## Find robust rare cell types
#' robust.clusters <- ConsensusRare(sc_object,
#'                                  assay = 'RNA',
#'                                  reduction = 'pca',
#'                                  dims = 1:50,
#'                                  k.param = 20,
#'                                  k = 6,
#'                                  Q_cut = 0.6,
#'                                  ratio = 0.2,
#'                                  reps = 30)
#'


ConsensusRare <- function(sc_object, assay='RNA', reduction='pca', dims=1:50, k.param=20, k=6, Q_cut=0.6, ratio=0.2, reps=30, user_num=NULL){

  # Step1: Run FindRare on shuffled data for multiple times and save results
  print('1. Start data shuffling:')
  n = dim(sc_object)[2]
  cluster_results <- matrix(NA, nrow = n, ncol = reps)

  cluster.num = c()
  for (r in 1:reps) {
	print(paste0('Data shuffling: ', r))
    set.seed(2^r)
    shuffled_idx <- sample(x = 1:n, size = n, replace = FALSE)

    sc_object_shuffle = sc_object
    sc_object_shuffle@reductions[[reduction]]@cell.embeddings = sc_object@reductions[[reduction]]@cell.embeddings[shuffled_idx,]
    sc_object_shuffle = FindNeighbors(sc_object_shuffle,
                                      reduction=reduction,
                                      k.param = k.param,
                                      compute.SNN = F,
                                      prune.SNN = 0,
                                      dims = dims,
                                      force.recalc = T,
                                      return.neighbor = T,
									  verbose = F)

    cluster_r <- FindRare(sc_object = sc_object_shuffle,
                                     assay = assay,
                                     k = k,
                                     Q_cut = Q_cut,
                                     ratio = ratio)

    cluster_results[, r] = cluster_r[order(shuffled_idx)]
    cluster_num = c(cluster.num, length(unique(cluster_r)))
  }

  # Step2: construct consensus matrix
  print("2. Build consensus matrix")
  consensus_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      both_present <- !is.na(cluster_results[i, ]) & !is.na(cluster_results[j, ])
      total <- sum(both_present)
      if (total == 0) {
        consensus_matrix[i, j] <- 0
      } else {
        consensus_matrix[i, j] <- sum(cluster_results[i, both_present] == cluster_results[j, both_present]) / total
      }
      consensus_matrix[j, i] <- consensus_matrix[i, j]
    }
  }
  diag(consensus_matrix) <- 1

  # Step3: clustering based on consensus matrix
  print('3. Output final cluster')
  freq <- table(cluster_num)
  cluster.num = as.integer(names(freq)[freq == max(freq)])
  print(paste0("The most frequent cluster numbers are: ", cluster.num))

  if(is.null(user_num)){
    k.num = cluster.num[1]
    print(paste0('The final cluster number is:', k.num))
  }else{
    k.num = user_num
    print('User defined cluster number is:', k.num)
  }

  final_cl <- hclust(as.dist(1 - consensus_matrix), method = "ward.D2")
  final_labels <- cutree(final_cl, k = k.num)

  return(final_labels)
}





