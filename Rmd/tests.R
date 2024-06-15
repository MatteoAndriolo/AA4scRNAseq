source("imports.R")
source("classes.R")

# Remove this part
TEST <- TRUE
HVF <- FALSE
TEST_genes <- 300
TEST_samples <- 500
CLASS.NAME <- "Exp1"
pathw <- "MAPK signaling pathway"
out_path <- sprintf("/app/out/%s/%s_files", CLASS.NAME, CLASS.NAME)

# Exp1
obj <- new("Exp1")
data_path <- "../data/AllonKleinLab/Experiment1/stateFate_inVitro_normed_counts.mtx"
cell_metadata <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_metadata.txt", header = TRUE, sep = "\t")
gene_names <- read.table("../data/AllonKleinLab/Experiment1/stateFate_inVitro_gene_names.txt")

se <- Matrix::readMM(data_path)
se@Dimnames[[2]] <- gene_names$V1
# se@Dimnames[[1]]=cell_metadata$Cell.barcode
# remove duplicated rows columns looking at names
# se <- se[!duplicated(se@Dimnames[[1]]), !duplicated(se@Dimnames[[2]])]
# se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
se <- CreateSeuratObject(counts = se)

if (!is.null(pathw)) {
  obj@pathw <- pathw
  obj <- obj_setGenes(obj, pathw)
  se <- se[gene_names, ]
} else if (test) {
  tgenes <- min(TEST_genes, nrow(se))
  tsamples <- min(TEST_samples, ncol(se))
  se <- se[1:tgenes, 1:tsamples]
  rm(tgenes, tsamples)
  se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
}

se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]
se <- ScaleData(se, layer = "counts")
se <- FindVariableFeatures(se)

se <- RunPCA(se, features = VariableFeatures(se))
se <- RunUMAP(se, features = VariableFeatures(se))

if (HVF) {
  m <- se@assays$RNA@layers$counts[which(se@assays$RNA@meta.data$vf_vst_counts_rank > 0), ]
} else {
  m <- se@assays$RNA@layers$counts
}

m <- m[Matrix::rowSums(m) > 0, Matrix::colSums(m) > 0]
m <- as.matrix(m)

obj@se <- se
obj@m <- m
