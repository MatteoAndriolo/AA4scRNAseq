source("imports.R")
source("classes.R")

# Remove this part
test <- FALSE
HVF <- FALSE
TEST_genes <- 300
TEST_samples <- 500
CLASS.NAME <- "Exp1"
pathw <- "MAPK signaling pathway"
out_path <- sprintf("/app/out/%s/%s_files", CLASS.NAME, CLASS.NAME)
obj <- new("Melanoma")

# read data
gene_names <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")
cell_metadata <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt", header = TRUE, sep = "\t")
new.names <- paste0(cell_metadata$Library, "_", cell_metadata$Cell.barcode)
cell_metadata$new.names <- new.names

if (!is.null(pathw)) {
  obj@pathw <- pathw
  obj <- obj_setGenes(obj, pathw)
  gene.flag <- gene_names$V1 %in% obj@genes
}

data_path <- "../data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx"
m <- Matrix::readMM(data_path)
m <- t(m)
m.org <- m

m <- m[gene.flag, ]
gene_names <- gene_names$V1[gene.flag]

if (is.null(pathw) & test) {
  tgenes <- max(TEST_genes, nrow(m))
  tsamples <- min(TEST_samples, ncol(m))
  tsamples <- 30000
}

se <- CreateSeuratObject(counts = m)
se <- ScaleData(se, layer = "counts")
se <- FindVariableFeatures(se)
se <- RunPCA(se, features = VariableFeatures(se))
se <- RunUMAP(se, features = VariableFeatures(se))

se.org <- se
dimnames(se) <- list(gene_names, new.names)
se@meta.data <- cell_metadata

if (HVF) {
  m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
} else {
  m <- se@assays$RNA@layers$counts
}

m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
m <- as.matrix(m)

obj@se <- se
obj@m <- m

obj_performArchetypes(obj)
