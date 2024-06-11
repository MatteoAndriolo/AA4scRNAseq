TH_archetypes.umap <- function(se, a, out_path=NULL,treshold=0.2){
  umap_result <- UMAPPlot(se)
  umap_data <- as.data.frame(umap_result$data)[,1:2]
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  
  column_sums <- colSums(a$archetypes)  # Sum of each column
  normalized_mat <- sweep(a$archetypes, 2, column_sums, FUN = "/")  # Divide each element by its column sum
  weights=as.data.frame(normalized_mat)
  
  plot_list <- list()
  for (i in 1:k) {
    umap_data$weight <- t(weights[i,])
    plot_title <- sprintf("UMAP Archetype %d", i)
    
    umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = weight)) +
      geom_point(size = 1) +
      scale_color_gradient(low = "grey", high = "red") +
      ggtitle(plot_title) +
      labs(color = "Weight")
    
    plot_list[[i]] <- umap_plot
    
    if(!is.null(out_path)){
      imgname <- sprintf("%s/UMAP_Archetype_%d.png", out_path, i)
      ggsave(imgname, plot = umap_plot)
      message(sprintf("Saving Image --- %s", imgname))
    }
  }
  
  # Combine all plots into a single image
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  
  # Save the combined image
  if(is.null(out_path)){
    combined_imgname <- sprintf("%s/UMAP_Combined.png", out_path)
    ggsave(combined_imgname, plot = combined_plot, width = 12, height = 8)
    message(sprintf("Saving Combined Image --- %s", combined_imgname))
  }
  
  return(combined_plot)
}