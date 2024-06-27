## AA Methods 
Test packages for performances and time
* [R-package](https://cran.r-project.org/web/packages/archetypes/index.html)
* [python-archetypes](https://github.com/aleixalcacer/archetypes)
* [python-archetypal](https://github.com/atmguille/archetypal-analysis)

## Clustering Methods
* Use k-means as reference
    * alternative Seurat????

## Datasets
|Type|Dim|Name|
|---|---|---|
|TIMESeries| piccolo | MouseCortex|
|TIMESeries |grande  | AllonKleinLab|
|Snapshot| piccolo| Melanoma|
|Snapshot|grande|Myocardite|


# Loading dataset flowchart
* If exists obj@se.org then copy it
* If it does not exists load data from 0
    * read data, gene_names and cell_metadata
    * (opt) remove row and columns with all 0 using `se <- se[Matrix::rowSums(se) > 0, Matrix::colSums(se) > 0]`
    * [EXP - DUPLICATES] remove duplicated genes cells
    * [MELANOMA - TEST] IF TEST
        * calculate tgenes, tsamples 
        * reduce se dimension
        * removing rows and columns all 0
    * CreateSeuratObject
    * save se.org
    * [EXP - DIMNAMES] 
        * give dimnames
        * AddMetadata
    * [EXP - TEST] IF TEST
        * calculate tgenes, tsamples
        * reduce dimension of obj@se
        * removing rows and columns all 0
    * Scale data 
    * FindVariableFeatures
    * RunPCA - features=rownames(obj@se)
    * RunUMAP - features=rownames(obj@se)
    * save obj@se.org
    * HVF - if HVF then reduce dimension of se to HVF anche check for 0 rows and columns

* if pathw
    * fetch pathw genes 
    * remove features not in pathw
    * !!!! Scale data 
    * !!!! FindVariableFeatures
    * !!!! RunPCA - features=rownames(obj@se)
    * !!!! RunUMAP - features=rownames(obj@se)


# Functions description

* obj_performArchetypes | obj | k [NULL] | doparallel [TRUE]
    * (param) if k is null calculate using elbowplot
    * 