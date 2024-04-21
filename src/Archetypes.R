library(archetypes)
library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(Single)

a=read.table("data/MelanomaSomethingPiccolo/GSE72056_melanoma_single_cell_revised_v2.txt")


meta_row_indices=2:5
cell_row_indices=5:nrow(a)
col_indices=2:ncol(a)
colNames=a[1,col_indices]

metadata=as.data.frame(a[2:4,col_indices], row.names = a[2:4,1])# Set the column names correctly
colnames(metadata) = colNames

rowNames=a[cell_row_indices,1]
ge=as.data.frame(a[cell_row_indices, col_indices], )


# using dplyr

# Create metadata data frame using dplyr
metadata <- a %>%
  slice(meta_row_indices) %>%
  select(col_indices) %>%
  `colnames<-`(colNames) %>%
  as.data.frame()


# ---------------------------------------------
metadata=as.matrix(a[2:4,2:ncol(a)])
rownames(metadata)=a[2:4,1]
colnames(metadata)=a[1,2:ncol(a)]

gene_expr=as.data.frame.matrix(a[5:nrow(a),2:ncol(a)], row.names=a[5:nrow(a),1])


gene_expr <- gene_expr %>% 
  mutate(across(where(is.character), as.numeric))
sp_gene_expr <- Matrix(gene_expr, sparse = TRUE)

# remove duplicate rows
rn=rownames(sp_gene_expr)
cn=unique(rownames(sp_gene_expr))
       
drn=rn[duplicated(rn)]
drn=which(rn==drns[1])
dcn=cn[duplicated(cn)]
gene_expr[drn[2],,drop=TRUE]


pbmc=CreateSeuratObject(sp_gene_expr,assay="scRNA")