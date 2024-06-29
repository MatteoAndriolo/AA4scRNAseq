out_folder="/app/out/Melanoma/Melanoma_files/uniHAnalysis5_12/"
save_slots_to_rds <- function(obj,out_folder) {
  # Get the class of the object
  obj_class <- class(obj)
  
  # Ensure the object is of class 'Melanoma'
  if (obj_class != "Melanoma") {
    stop("The object is not of class 'Melanoma'")
  }
  
  # Get the names of all slots
  slot_names <- slotNames(obj)
  
  # Iterate through each slot
  for (slot_name in slot_names) {
    # Extract the content of the slot
    slot_content <- slot(obj, slot_name)
    
    # Create the file name
    file_name <- paste0(out_folder,slot_name, ".rds")
    
    # Save the slot content to the file
    saveRDS(slot_content, file_name)
    
    cat("Saved slot", slot_name, "to", file_name, "\n")
  }
}

# Define the function to save each plot in the 'plots' slot to a separate .png file
save_plots_to_png <- function(obj, output_folder) {
  # Get the class of the object
  obj_class <- class(obj)
  
  # Ensure the object is of class 'Melanoma'
  if (obj_class != "Melanoma") {
    stop("The object is not of class 'Melanoma'")
  }
  
  # Check if the 'plots' slot exists
  if (!("plots" %in% slotNames(obj))) {
    stop("The object does not have a 'plots' slot")
  }
  
  # Extract the plots
  plots <- slot(obj, "plots")
  
  # Ensure the 'plots' slot contains multiple plots
  if (!is.list(plots)) {
    stop("The 'plots' slot does not contain a list of plots")
  }
  
  # Create the output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Iterate through each plot and save as a PNG
  for (i in seq_along(plots)) {
    plot <- plots[[i]]
    # Create the file name
    file_name <- file.path(output_folder, paste0("plot_", i, ".png"))
    
    # Save the plot as a PNG
    png(file_name)
    print(plot)
    dev.off()
    
    cat("Saved plot", i, "to", file_name, "\n")
  }
}


source("/app/Rmd/imports.R")
source("/app/Rmd/classes.R")
obj <- readRDS("/app/out/Melanoma/Melanoma_files/uniHAnalysis5_12/MelanomaH_5.13.FS.unique..rds")
#save_slots_to_rds(obj, out_folder)
save_plots_to_png(obj, out_folder)
