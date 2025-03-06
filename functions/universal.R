#Plot MADs for diff data in a seurat object---------
outlierHistogramSeurat <- function(sobj, feature, mads = 1:3, bins = 30,
                                   show_zero = FALSE) {
  # Extract metadata from the Seurat object
  metadata <- sobj@meta.data

  # Check if the feature exists in metadata
  if (!(feature %in% colnames(metadata))) {
    stop(paste("Feature", feature, "not found in the Seurat object's metadata"))
  }

  # Extract the feature data from the metadata
  feature_data <- metadata[[feature]]

  # Check if the feature data is numeric
  if (!is.numeric(feature_data)) {
    stop(paste("Feature", feature, "is not numeric. Please check your data."))
  }

  # Calculate median and MAD (Median Absolute Deviation)
  med <- median(feature_data)
  MAD <- mad(feature_data, center = med, na.rm = TRUE)

  # Identify outliers using the MAD method
  outs <- lapply(mads, function(n) {
    lower <- med - n * MAD
    higher <- med + n * MAD
    n_low  <- sum(feature_data < lower)
    n_high <- sum(feature_data > higher)

    list(n = n, lower = lower, higher = higher,
         n_low = n_low, n_high = n_high)
  })

  # Plot the histogram of the feature
  gg <- ggplot(data.frame(feature_data), aes(x = feature_data)) +
    geom_histogram(bins = bins, fill = "grey", color = "black") +
    geom_vline(xintercept = med, colour = "blue") +
    annotate("text", x = med, y = 0, label = "Median", colour = "blue",
             angle = 90, hjust = -0.1, vjust = -0.5) +
    annotate("text", x = med, y = Inf, label = round(med, 2), colour = "blue",
             angle = 90, hjust = 1, vjust = -0.5) +
    theme_minimal() +
    theme(legend.position = "none")

  # Define a color palette for MAD lines
  cols <- scales::brewer_pal(palette = "Reds", direction = -1)(length(mads) + 1)

  # Add lines for the outlier thresholds
  for (i in seq_along(outs)) {
    out <- outs[[i]]

    if (show_zero) {
      show_low  <- TRUE
      show_high <- TRUE
    } else {
      show_low  <- out$n_low > 0
      show_high <- out$n_high > 0
    }

    if (show_low) {
      gg <- gg +
        geom_vline(xintercept = out$lower, colour = cols[i]) +
        annotate("text", x = out$lower, y = 0,
                 label = paste0("-", out$n, " MADs"), colour = cols[i],
                 angle = 90, hjust = -0.1, vjust = -0.5) +
        annotate("text", x = out$lower, y = Inf,
                 label = paste0(round(out$lower, 2),
                                " (", out$n_low, " lower)"),
                 colour = cols[i], angle = 90, hjust = 1, vjust = -0.5)
    }

    if (show_high) {
      gg <- gg +
        geom_vline(xintercept = out$higher, colour = cols[i]) +
        annotate("text", x = out$higher, y = 0,
                 label = paste0("+", out$n, " MADs"), colour = cols[i],
                 angle = 90, hjust = -0.1, vjust = -0.5) +
        annotate("text", x = out$higher, y = Inf,
                 label = paste0(round(out$higher, 2),
                                " (", out$n_high, " higher)"),
                 colour = cols[i], angle = 90, hjust = 1, vjust = -0.5)
    }
  }

  return(gg)
}

#example
#outlierHistogramSeurat(sobj, feature="nCount_RNA", mads = 1:6)

#Extract all saved files from a script-----------
extract_saved_files <- function(script_path) {
  # Read the script file
  script_lines <- readLines(script_path)

  # Regular expressions to match save() and saveRDS() calls
  save_pattern <- "save\\((.*?),(.*?)\\)"  # Matches save(file, ...)
  save_rds_pattern <- "saveRDS\\((.*?),\\s*\"(.*?)\"\\)"  # Matches saveRDS(obj, "file")

  # Initialize an empty data frame
  results <- data.frame(
    File = character(),
    Directory = character(),
    Description = character(),
    stringsAsFactors = FALSE
  )

  # Extract save() and saveRDS() calls
  for (line in script_lines) {
    save_match <- regmatches(line, regexec(save_pattern, line))[[1]]
    save_rds_match <- regmatches(line, regexec(save_rds_pattern, line))[[1]]

    # If save() pattern is found
    if (length(save_match) > 0) {
      file_name <- gsub("\"|'", "", save_match[2])  # Extract the file name
      directory <- gsub("\"|'", "", save_match[3])  # Extract the directory
      description <- readline(prompt = paste("Enter a description for the directory:", directory, ": "))
      results <- rbind(results, data.frame(
        File = file_name,
        Directory = directory,
        Description = description,
        stringsAsFactors = FALSE
      ))
    }

    # If saveRDS() pattern is found
    if (length(save_rds_match) > 0) {
      file_name <- gsub("\"|'", "", save_rds_match[2])  # Extract the file name
      directory <- file_name  # For saveRDS, the directory is just the file path
      description <- readline(prompt = paste("Enter a description for the directory:", directory, ": "))
      results <- rbind(results, data.frame(
        File = file_name,
        Directory = directory,
        Description = description,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(results)
}

#example
#saved_files <- extract_saved_files("./analysis/02_quality_control.qmd")
#> head(saved_files)
# File                                                                      Directory
# 3      selected                          file=here::here(output/processed/02_SI_selected.Robj)
# 4          sobj        file=/home/nsingh/sterile_inflammation/output/processed/02_SI_sobj.Robj
# 5          sobj  file=/home/nsingh/sterile_inflammation/output/processed/02_SI_sobj_scale.Robj
# 6 filtered_sobj                                     file=./output/processed/02_SI_sobj_qc.Robj
# 7 filtered_sobj     file=/home/nsingh/sterile_inflammation/output/processed/02_SI_sobj_df.Robj
# 8  sobj_post_df                                file=./output/processed/02_SI_sobj_post_df.Robj

#Function to calculate how many PCs to use for UMAP-----------
select_dims <- function(object, reduction) {
  # Approach 1:
  # The point where individual principal components only contribute 5% of standard deviation
  # and the principal components cumulatively contribute 90% of the standard deviation.
  pct <- object[[reduction]]@stdev / sum(object[[reduction]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  print(paste('Approach 1: Dims = ', co1, sep = ''))

  # Approach 2:
  # The point where the percent change in variation between consecutive PCs is < 0.1%.
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  print(paste('Approach 2: Dims = ', co2, sep = ''))

  # Approach 3:
  # The point where the cumulative variance explained by the PCs is at least 95%.
  # We square the PC standard deviations to get variance, sum them, and compute cumulative fraction.
  stdevs <- object[[reduction]]@stdev
  var_per_pc <- stdevs^2
  total_variance <- sum(var_per_pc)
  var_explained <- var_per_pc / total_variance
  cum_var_explained <- cumsum(var_explained)
  co3 <- which(cum_var_explained >= 0.95)[1]
  print(paste('Approach 3: Dims = ', co3, sep = ''))

  # Print a brief explanation of each approach
  cat("\nSummary of Approaches:\n")
  cat("Approach 1: The point where PCs collectively contribute ~90% of the standard deviation and the individual PCs only contribute 5% of standard deviation.\n")
  cat("Approach 2: The point where the percent change in variation between consecutive PCs is < 0.1%.\n")
  cat("Approach 3: The point where the cumulative variance explained by the PCs is at least 95%.\n\n")

  # Return the dims based on the minimum of the three approaches
  return(1:min(co1, co2, co3))
}

# Define the function to add mini arrow axes-------------
# function to add_mini_axis() that you can use in conjunction with your DimPlot (or any ggplot object) to overlay mini arrow axes.
# Load necessary libraries (if not already loaded)
library(ggplot2)
library(cowplot)
library(grid)

add_mini_axis <- function(p,
                          x = 0.01, y = 0.01,
                          width = 1, height = 1,
                          arrow_length = unit(0.1, "inches"),
                          arrow_type = "open",
                          arrow_angle = 25,
                          label_fontsize = 10) {

  # Define arrow style using grid's arrow function
  arrow_style <- arrow(length = arrow_length, type = arrow_type, angle = arrow_angle)

  # Create a custom grob for the mini axes with arrows
  axis_grob <- grobTree(
    # Horizontal arrow for UMAP1
    segmentsGrob(
      x0 = unit(0.05, "npc"), x1 = unit(0.25, "npc"),
      y0 = unit(0.05, "npc"), y1 = unit(0.05, "npc"),
      gp = gpar(col = "black", lwd = 2),
      arrow = arrow_style
    ),
    # Vertical arrow for UMAP2
    segmentsGrob(
      x0 = unit(0.05, "npc"), x1 = unit(0.05, "npc"),
      y0 = unit(0.05, "npc"), y1 = unit(0.25, "npc"),
      gp = gpar(col = "black", lwd = 2),
      arrow = arrow_style
    ),
    # Label for UMAP1
    textGrob("UMAP1", x = unit(0.15, "npc"), y = unit(0.03, "npc"),
             gp = gpar(fontsize = label_fontsize)),
    # Label for UMAP2 (rotated so it appears vertically)
    textGrob("UMAP2", x = unit(0.03, "npc"), y = unit(0.15, "npc"),
             rot = 90, gp = gpar(fontsize = label_fontsize))
  )

  # Overlay the mini axes on the provided ggplot object using cowplot's ggdraw
  p_final <- ggdraw(p) +
    draw_grob(axis_grob, x = x, y = y, width = width, height = height)

  return(p_final)
}

# Example usage:
# Suppose you have a DimPlot stored in p4:
# p4 <- DimPlot(sobj, reduction = "umap") + theme_minimal()
# Now add the mini axes:
# p_final <- add_mini_axis(p4)
# print(p_final)

#Function to add title to featureplots---------------
library(Seurat)
library(cowplot)

make_feature_panel <- function(seurat_obj, gene_list, panel_title, ncol = 2) {
  # Create individual FeaturePlots for each gene without combining them
  p_list <- FeaturePlot(seurat_obj, features = gene_list, combine = FALSE,
                        label.size = 1)

  #Remove axes
  p_list <- lapply(p_list, function(p) {
    p + theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
  })

  # Arrange the individual plots into a grid with the specified number of columns
  grid <- plot_grid(plotlist = p_list, ncol = ncol) + theme_void()

  # Add an overall panel title above the grid
  final_plot <- ggdraw() +
    draw_label(panel_title, fontface = "bold", size = 16, x = 0.5, y = 0.98, hjust = 0.5) +
    draw_plot(grid, y = 0, height = 0.95)

  return(final_plot)
}

# Example usage:
# # Specify your gene list and panel title
# endothelial_genes <- c("Pecam1", "Cd34", "Ace", "Vcam1")
# panel_title <- "Endothelial Cells"
#
# # Create the combined feature plot panel
# p_endothelial <- make_feature_panel(seurat_obj, endothelial_genes, panel_title, ncol = 2)
#
# # Display the plot
# print(p_endothelial)


