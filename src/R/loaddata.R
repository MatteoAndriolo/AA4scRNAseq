require(Matrix)
require(dplyr)
require(Seurat)
# require(monocle3)

TEST=FALSE
TEST_genes=3000
TEST_samples=5000

# SMALL SNAPSHOT
loadMelanoma <- function(hvf.nfeatures = 2000 ){
	# -------- BASIC SETUP
	ge <- read.table("data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt", header = TRUE)
	ge=ge[!duplicated(ge[,1]),] # remove duplicated genes
	rownames(ge)=ge[,1]         # extract from matrix rownames
	ge=ge[,2:ncol(ge)]          # elide rownames from gene expression matrix

	metadata=ge[1:3,]          # extract metadata
	metadata <- t(metadata) %>%
		data.frame() %>% 
		mutate(across(where(is.character), as.numeric)) 

	ge=ge[4:nrow(ge),]
	ge <- ge %>% 
		data.frame() %>% 
		mutate(across(where(is.character), as.numeric)) 
	# -------- END BASIC SETUP

	if(TEST){
		tgenes=min(TEST_genes, nrow(ge))
		tsamples=min(TEST_samples, ncol(ge))
		metadata=metadata[1:tsamples,]
		ge=ge[1:tgenes,1:tsamples]
		rm(tgenes,tsamples)
	}

	# -------- SEURAT
	se <- CreateSeuratObject(ge,project="MelanomaSC",meta.data=metadata)
	# Skip normalization since data is already normalized
	se <- ScaleData(se, layer= "counts") #,check.for.norm=FALSE)
	se <- FindVariableFeatures(se, hvf.nfeatures = hvf.nfeatures)

	se = RunPCA(se, features = VariableFeatures(se))
	se = RunUMAP(se, features = VariableFeatures(se))

	return(se)
}

# LARGE SNAPSHOT
loadMyocardialInfarction <- function(hvf.nfeatures = 2000){
	se <- readRDS("data/MyocardialInfarction/e61af320-303a-4029-8500-db6636bba0d4.rds")
	se <- FindVariableFeatures(se)
	warning("Implementation Not completed due to problem with data")
	stop()

	se = RunPCA(se, features = VariableFeatures(se))
	se = RunUMAP(se, features = VariableFeatures(se))

	return(se)
}

# SMALL TIMESERIES
loadMouseCortex <- function(hvf.nfeatures=2000){
	load("data/MouseCortex/MouseCortex.RData")
	MouseCortex=MouseCortex
	warning("Implementation Not completed due to problem with data")

	se = RunPCA(se, features = VariableFeatures(se))
	se = RunUMAP(se, features = VariableFeatures(se))

	stop()
}

# LARGE TIMESERIES
loadExp1 <- function(hvf.nfeatures = 2000){
	#_ # Binary matrix indicating clonal membership of each cell
	#_ # The rows of this file represent cells and correspond to the rows of _counts_matrix_in_vitro_ (above).
	#_ # The columns represent clones. Not every cell belongs to a clone. 
	#_ clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment1/stateFate_inVitro_clone_matrix.mtx")
	#_  # cell metadata : cell type annotation
	#_  cell_metadata <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt",header=TRUE,sep = "\t" )
	#_  
	#_  gene_names <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")
	#_  
	#_  # List of cells belonging to the neutrophil/monocyte trajectory that were used in becnmark analysis
	#_  neutrophil_monocyte_trajectory <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_monocyte_trajectory.txt",header=TRUE,sep="\t")
	#_  # pseudotime for neutrophil trajectory cells
	#_  neutrophil_pseudotime <- read.table("data/AllonKleinLab/Experiment1/stateFate_inVitro_neutrophil_pseudotime.txt",header=TRUE, sep="\t")

	#  Number of transcripts for each gene in each cell
	se <- Matrix::readMM("/app/data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx" )

	if(TEST){
		tgenes=min(TEST_genes, nrow(se))
		tsamples=min(TEST_samples, ncol(se))
		se=se[1:tgenes,1:tsamples]
		rm(tgenes,tsamples)
	}

	se <- CreateSeuratObject(counts=se)
	se <- ScaleData(se, layer= "counts") #,check.for.norm=FALSE)
	se <- FindVariableFeatures(se)

	se = RunPCA(se, features = VariableFeatures(se))
	se = RunUMAP(se, features = VariableFeatures(se))
	
	return(se)
}

loadExp2 <- function(hvf.nfeatures=2000){
	# clone_matrix <- Matrix::readMM("data/AllonKleinLab/Experiment2/stateFate_inVivo_clone_matrix.mtx")
	# gene_names <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_gene_names.txt",sep="\t")
	# metadata <- read.table("data/AllonKleinLab/Experiment2/stateFate_inVivo_metadata.txt",sep="\t")

	se <- Matrix::readMM("data/AllonKleinLab/Experiment2/stateFate_inVivo_normed_counts.mtx")

	if(TEST){
		tgenes=min(TEST_genes, nrow(se))
		tsamples=min(TEST_samples, ncol(se))
		se=se[1:tgenes,1:tsamples]
		rm(tgenes,tsamples)
	}

	se <- CreateSeuratObject(counts=se)
	se <- ScaleData(se, layer="counts")
	se <- FindVariableFeatures(se)

	se = RunPCA(se, features = VariableFeatures(se))
	se = RunUMAP(se, features = VariableFeatures(se))

	return(se) 
}

loadExp3 <- function(hvf.nfeatures = 2000){
	# clone_matrix = Matrix::readMM("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_clone_matrix.mtx")
	# gene_names <- read.table( "data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_gene_names.txt")
	# metadata <- read.table("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_metadata.txt", sep="\t")

	se <- Matrix::readMM("data/AllonKleinLab/Experiment3/stateFate_cytokinePerturbation_normed_counts.mtx")

	if(TEST){
		tgenes=min(TEST_genes, nrow(se))
		tsamples=min(TEST_samples, ncol(se))
		se=se[1:tgenes,1:tsamples]
		rm(tgenes,tsamples)
	}

	se <- CreateSeuratObject(counts=se)
	se <- ScaleData(se, layer= "counts") #,check.for.norm=FALSE)
	se <- FindVariableFeatures(se)

	se = RunPCA(se, features = VariableFeatures(se))
	se = RunUMAP(se, features = VariableFeatures(se))

	return(se)
}
