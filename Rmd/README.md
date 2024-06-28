# parameter for unique.Rmd <- _main.sh <- factory???

```R
# Define the 'database' Class
setClass("database",
  slots = list(
    se = "Seurat",                 # Seurat object
    a = "Archetypes",              # Archetype model (assuming a specific class for archetypes)
    m = "matrix",                  # Data matrix
    umap.archetypes = "ggplot",    # UMAP plot for archetypes
    combined_plot = "ggplot",      # Combined PCA and UMAP plots
    elbowplot = "ggplot",          # Elbow plot for PCA
    umap = "ggplot",               # UMAP plot
    pca = "ggplot",                # PCA plot
    simplexplot = "ggplot",        # Plot for simplex representation of archetypes
    pathw = "character",           # Pathway information
    genes = "character"            # List of genes
  )
)
```