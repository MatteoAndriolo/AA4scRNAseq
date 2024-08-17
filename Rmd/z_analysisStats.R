source("Rmd/imports.R")
source("Rmd/classes.R")
source("Rmd/class_Melanoma.R")
source("Rmd/class_Mouse.R")
set.seed(2024)
NOT_FINAL <- TRUE
debug <- TRUE

obj <- new("Melanoma")
obj <- new("Mouse")


# for (pw in list("FS1")) {
for (pw in list("FS1", "FS2", "FS3", "FS4", "FS5", "HFS")) {
  if (class(obj) == "Melanoma") {
    obj@params$out_path <- paste("outPar/Melanoma/0813_0758/", pw, "_1738961", sep = "")
  } else if (class(obj) == "Mouse") {
    obj@params$out_path <- paste("outPar/Mouse/0813_0758/", pw, "_1738962", sep = "")
  }
  if (pw == "HFS") {
    obj@params$hvf <- TRUE
    obj@params$pathw <- NULL
    obj@other$namePathw <- "HVF"
  } else {
    obj@params$hvf <- FALSE
    obj@params$pathw <- as.integer(substr(pw, 3, 3))
    obj@other$namePathw <- list(
      "GLYK",
      "MAPK",
      "CANCER",
      "MTOR",
      "TGF"
    )[[obj@params$pathw]]
  }

  obj@params$path_outdata <- file.path(obj@params$out_path, "data")

  # insefile <- file.path(obj@params$path_outdata, paste(ifelse(pw == "HFS", "", substr(pw, 3, 3)), ".Rds", sep = ""))
  # obj@se <- readRDS(insefile)

  inaafile <- file.path(obj@params$path_outdata, "archetypes.Rds")
  obj@archetypes <- readRDS(inaafile)

  inmdfile <- file.path(obj@params$path_outdata, "metadata.Rds")
  obj@other <- readRDS(inmdfile)

  message(ifelse(obj@other$namePathw == "HVF", "HVF", paste(obj@params$pathw, obj@other$namePathw, sep = "")))
  obj@params$path_figures <- file.path(
    obj@params$out_path,
    paste(
      ifelse(class(obj) == "Melanoma", "MEL", "MOUSE"),
      "_",
      ifelse(obj@other$namePathw == "HVF", "HVF", paste(obj@params$pathw, obj@other$namePathw, sep = "")),
      sep = ""
    )
  )
  obj@other$treshold <- 0.5



  ##################################################
  # BEGIN WITH FIGURES
  ##################################################
  plot_data <- bind_rows(obj@archetypes$analysis, .id = "narch") %>%
    mutate(narch = as.integer(narch))
  head(plot_data)
  #                sse        varexpt   time num_archetypes
  # 6  463667.438314 0.758519932482  40.38              6
  # 7  452990.382946 0.764080590484  65.36              7
  # 8  440877.679379 0.770388940464 102.59              8
  # 9  430837.696003 0.775617808534  34.49              9
  # 10 422259.030129 0.780085615939  28.20             10
  # 11 413308.365342 0.784747162036  32.80             11

  # Plot times
  # use plot_data$num_archetypes as x-axis
  p <- ggplot(plot_data, aes(x = narch, y = time)) +
    geom_point() +
    labs(
      x = "Number of archetypes",
      y = "Time (s)"
    ) +
    scale_x_continuous(breaks = min(plot_data$narch):max(plot_data$narch))
  theme_classic()

  ggsave(file.path(obj@params$path_figures, "AA_time.png"), p, width = 8, height = 6)

  # Plot SSE
  p <- ggplot(plot_data, aes(x = narch, y = sse)) +
    geom_point() +
    labs(
      x = "Number of archetypes",
      y = "SSE"
    ) +
    scale_x_continuous(breaks = min(plot_data$narch):max(plot_data$narch))
  theme_classic()

  ggsave(file.path(obj@params$path_figures, "AA_sse.png"), p, width = 8, height = 6)

  # Plot Varexpt
  p <- ggplot(plot_data, aes(x = narch, y = varexpt)) +
    geom_point() +
    labs(
      x = "Number of archetypes",
      y = "Varexpt"
    ) +
    scale_x_continuous(breaks = min(plot_data$narch):max(plot_data$narch))

  ggsave(file.path(obj@params$path_figures, "AA_varexpt.png"), p, width = 8, height = 6)
}

##################################################
# SANKEY PLOT
##################################################
for (treshold in c(TRUE, FALSE)) {
  # Create a data frame with the required columns
  View(obj@other$se.metadata)
  if (treshold) {
    data <- as.data.frame(list(
          type = plot_data$ctype,
          archetype = plot_data$aaclusters.treshold
        ))
  } else {
    data <- as.data.frame(list(
          type = plot_data$ctype,
          archetype = plot_data$aaclusters
        ))
  }

  data$archetype <- factor(
        data$archetype,
        levels = c("1", "2", "3", "4", "5", "6", "7", "NotAssigned", "Archetype"),
        labels = c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "NotAssigned", "Archetype")
      )

  d <- data %>%
        make_long(colnames(data)) %>%
        filter(node != "Archetype") # %>% filter(next_node != "Archetype")


  colorMapArchetypesSankey <- setNames(
        colorMapArchetypes,
        c(paste0("A", 1:num_archetypes, sep = ""), "NotAssigned", "Archetype")
      )

  # END SETUP #################################
  # BASE SANKEY
  if (TRUE) {
    pl <- ggplot(d, aes(
          x = x,
          next_x = next_x,
          node = node,
          next_node = next_node,
          fill = factor(node),
          label = node
        )) +
          geom_sankey(
            flow.alpha = 0.5,
            node.color = "black",
            show.legend = FALSE
          ) +
          geom_sankey_label(size = 3, color = "black", fill = "white") +
          scale_fill_manual(
            values = c(colorMapCTypes, colorMapArchetypesSankey)
          ) +
          scale_x_discrete(
            labels = c("type" = "Cell Type", "archetype" = "Archetype")
          ) +
          theme_alluvial() +
          theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank()
          )
    pl

    prefixName <- paste(obj@other$namePathw, k, ifelse(treshold > 0, "th", ""), sep = ".")
    ggsave(
          file.path(obj@params$path_figures, paste(prefixName, "sankey", ifelse(treshold > 0, "th", ""), "png", sep = ".")),
          pl,
          width = 8,
          height = 6
        )

    if (class(obj) == "Melanoma") {
      for (mal in c(0, 1, 2)) {
        which.is.mal <- which(plot_data$malignant == mal)
        if (treshold) {
          data <- as.data.frame(list(
                type = plot_data$ctype[which.is.mal],
                archetype = plot_data$aaclusters.treshold[which.is.mal]
              ))
        } else {
          data <- as.data.frame(list(
                type = plot_data$ctype[which.is.mal],
                archetype = plot_data$aaclusters[which.is.mal]
              ))
        }

        data$archetype <- factor(
              data$archetype[which.is.mal],
              levels = c("1", "2", "3", "4", "5", "6", "7", "NotAssigned", "Archetype"),
              labels = c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "NotAssigned", "Archetype")
            )

        d <- data %>%
              make_long(colnames(data)) %>%
              filter(node != "Archetype") # %>% filter(next_node != "Archetype")


        pl <- ggplot(d, aes(
              x = x,
              next_x = next_x,
              node = node,
              next_node = next_node,
              fill = factor(node),
              label = node
            )) +
              geom_sankey(
                flow.alpha = 0.5,
                node.color = "black",
                show.legend = FALSE
              ) +
              geom_sankey_label(size = 3, color = "black", fill = "white") +
              scale_fill_manual(
                values = c(colorMapCTypes, colorMapArchetypesSankey)
              ) +
              scale_x_discrete(
                labels = c("type" = "Cell Type", "archetype" = "Archetype")
              ) +
              theme_alluvial() +
              theme(
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line.x = element_blank(),
                axis.ticks.x = element_blank()
              )
        pl

        namesMalignant = c("Unknown", "Non.Malignant", "Malignant")
        prefixName <- paste(obj@other$namePathw, k, ifelse(treshold > 0, "th", ""), sep = ".")
        ggsave(
              file.path(obj@params$path_figures, paste(prefixName, "sankey", namesMalignant, "png", sep = ".")),
              pl,
              width = 8,
              height = 6
            )
      }
    }
  }

  # HEATMAP
  # remove from table all datapoints with Archetype and also factors (i dont wont row and columns to 0)
  # df <- table(data$type[-"Archetype"], data$archetype[-"Archetype"])
  df <- table(data$type, data$archetype)
  df <- df[, colnames(df) != "Archetype"]
  df <- df[row.names(df) != "Archetype",]
  df <- as.data.frame(df)

  # Heatmap with text inside
  plt_hm <- ggplot(df, aes(x = Var1, y = Var2)) +
        geom_tile(aes(fill = Freq), color = "white") +
        geom_text(aes(label = Freq), vjust = 1) +
        scale_fill_gradient(low = "white", high = "blue") +
        theme_alluvial() +
        labs(x = "Cell types", y = "Archetype") +
        theme(axis.text.x = element_text(hjust = 1))
  plt_hm

  ggsave(
        file.path(obj@params$path_figures, paste(prefixName, "heatmap", ifelse(treshold > 0, "th", ""), "png", sep = ".")),
        width = 8,
        height = 6
      )
}
# end treshold 2
