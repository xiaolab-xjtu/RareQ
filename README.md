# RareQ
RareQ (R package) is a network-propagation-based algorithm for rare cell type identification, that demonstrates superior performance and higher scalability to existing algorithms.


## Required R modules
```R
    R >= 4.0.0
    Seurat >= 4.0.2
    Signac >= 1.9.0  # For preprocessing scATAC-seq data
```

## Installation
```R
  # Install in R with devtools
  library(devtools)
  install_github('fabotao/RareQ')
```

## Usage
Example [Jurkat](https://github.com/fabotao/RareQ/blob/main/data/Jurkat.RDS) scRNA-seq data for demonstration

```R
  library(RareQ)
  library(Seurat) 
  
  # Read example data
  obj = readRDS('Jurkat.RDS)  # Example data from data folder
  counts = obj@assays$RNA@counts
  
  # Preprocessing scRNA-seq data
  sc_object <- CreateSeuratObject(count=counts, project = "sc_object", min.cells = 3)
  sc_object <- NormalizeData(sc_object)
  sc_object <- FindVariableFeatures(sc_object, nfeatures=2000)
  sc_object <- ScaleData(sc_object)
  sc_object <- RunPCA(sc_object, npcs=50)
  sc_object <- RunUMAP(sc_object, dims=1:50)
  sc_object <- FindNeighbors(object = sc_object,
                             k.param = 20, 
                             compute.SNN = F, 
                             prune.SNN = 0, 
                             reduction = "pca", 
                             dims = 1:50, 
                             force.recalc = F, 
                             return.neighbor = T)
  
  # Use RareQ to derive both major and rare cell clusters  
  cluster <- FindRare(sc_object)
  table(cluster)
  
  sc_object$cluster = cluster
  DimPlot(sc_object, group.by='cluster')
  
```

### (Optional)
Though deterministic label propagation in **RareQ** has demonstrated high robustness to data shuffling, we also provide a helper function **ConsensusRare** which runs FindRare on multiple shuffled datasets to derived more robust result via consensus clustering. 
**Attention:** It takes more time than FindRare function and is less memory efficient for large datasets due to confusion matrix.

```R
  library(RareQ)
  library(Seurat) 
  
  # Read example data
  obj = readRDS('Jurkat.RDS)  # Example data from data folder
  counts = obj@assays$RNA@counts
  
  # Preprocessing scRNA-seq data
  sc_object <- CreateSeuratObject(count=counts, project = "sc_object", min.cells = 3)
  sc_object <- NormalizeData(sc_object)
  sc_object <- FindVariableFeatures(sc_object, nfeatures=2000)
  sc_object <- ScaleData(sc_object)
  sc_object <- RunPCA(sc_object, npcs=50)
  sc_object <- RunUMAP(sc_object, dims=1:50)
  
  # Use ConsensusRare to derive more robust result  
  cluster <- ConsensusRare(sc_object,
                                  assay = 'RNA',      # Use RNA assay for default  
                                  reduction = 'pca',  # Use pca reduction
                                  dims = 1:50,
                                  k.param = 20,       # k.param typical in scRNA workflows
                                  k = 6,              # k parameter for FindRare
                                  Q_cut = 0.6,        # Q_cut parameter for FindRare                                 
								  ratio = 0.2,        # ratio paramter for FindRare
                                  reps = 30)          # Number of data reshuffling
)
  table(cluster)
  
  sc_object$cluster = cluster
  DimPlot(sc_object, group.by='cluster')
```

## Tutorial
We have tutorials to assist you in utilizing RareQ. You can locate these tutorials in the RareQ/Tutorial_example directory of the project. Additionally, we have provided related Dataset to aid you in testing and acquainting yourself with the RareQ functionality:

[Example Datasets](https://zenodo.org/records/17190972/files/Tutorial_example.rar?download=1)

a. The tutorial contains four R scripts encompassing different data modality as follows:
1. [scRNA_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/scRNA_analysis.html): RareQ analysis of scRNA-seq data from Jurkat cell line
2. [scRNA_scATAC_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/scRNA_scATAC_analysis.html): RareQ analysis of scRNA-seq + scATAC-seq multiome data from mouse T cells
3. [scRNA_ADT_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/scRNA_ADT_analysis.html): RareQ analysis of CITE-seq data (containing both RNA and ADT modalities) from human bone marrow mononuclear cells.
4. [Xenium_spatial_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/Xenium_spatial_analysis.html): RareQ analysis of Xenium spatial data from 10X Genomics from mouse brain

b. (Optional) The tutorial for ConsensusRare is as follows:
1. [scRNA_analysis](https://xiaolab-xjtu.github.io/RareQ/Tutorials/scRNA_analysis.html): ConsensusRare analysis of scRNA-seq data from Jurkat cell line

## Simulation
We have provided code and data to assist you in generating simulation datasets to benchmark RareQ. You can locate the Python script in the RareQ/Simulation directory of the project. 

[Datasets for simulation](https://zenodo.org/records/17190972/files/simulation.rar?download=1)


## Citation


