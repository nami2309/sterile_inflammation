#' Outlier histogram
#'
#' Histogram with lines showing median and MADs from median
#'
#' @param seurat_obj seurat_obj to plot
#' @param x Name of the metadata column to plot
#' @param mads Vector of mads to show
#' @return ggplot object
outlierHistogramSeurat <- function(seurat_obj, feature, mads = 1:3, bins = 30,
                                   show_zero = FALSE) {
  # Extract metadata from the Seurat object
  metadata <- seurat_obj@meta.data

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


# Example usage with Seurat object:
# outlierHistogramSeurat(seurat_obj, feature = "nCount_RNA")
# outlierHistogramSeurat(seurat_obj, feature = "percent.mt")
# outlierHistogramSeurat(seurat_obj, feature = "nCount_RNA", mads = 1:3)


#' PAGA cluster graph
#'
#' Plot PAGA cluster results
#'
#' @param embedding data.frame with positions for nodes
#' @param edges data.frame with edges
#' @param thresh threshold for connectivity
#' @param colour column to use for node colour
#'
#' @return ggplot object
plotPAGAClustGraph <- function(embedding, edges, thresh = 0,
                               colour = "Cluster") {

    is_discrete <- is.factor(embedding[[colour]])

    gg <- ggplot(embedding, aes(x = X, y = Y))

    if (is_discrete) {
        gg <- gg +
            geom_segment(data = filter(edges, Connectivity > thresh),
                         aes(x = FromX, y = FromY, xend = ToX, yend = ToY,
                             colour = Connectivity),
                         size = 4) +
            scale_colour_viridis_c(limits = c(0, 1))
    } else {
        gg <- gg +
            geom_segment(data = filter(edges, Connectivity > thresh),
                         aes(x = FromX, y = FromY, xend = ToX, yend = ToY,
                             alpha = Connectivity),
                         size = 4, colour = "grey30") +
            scale_alpha(limits = c(0, 1)) +
            scale_fill_viridis_c()
    }

    gg <- gg +
        geom_point(aes(fill = !!ensym(colour), size = Size), shape = 21) +
        geom_text(aes(label = Cluster)) +
        scale_size(range = c(5, 15)) +
        theme_void() +
        theme(legend.position = "none")

    return(gg)
}


#' PAGA cell graph
#'
#' Plot PAGA cell results
#'
#' @param embedding data.frame with positions for nodes
#' @param edges data.frame with edges
#' @param thresh threshold for connectivity
#' @param colour column to use for cell colour
#' @param label whether to add cluster labels
#'
#' @return ggplot object
plotPAGACellGraph <- function(embedding, edges, thresh = 0, colour = "Cluster",
                              label = TRUE) {

    is_discrete <- is.factor(embedding[[colour]])

    gg <- ggplot(embedding, aes(x = X, y = Y, colour = !!ensym(colour))) +
        geom_segment(data = filter(edges, Connectivity > thresh),
                     aes(x = FromX, y = FromY, xend = ToX, yend = ToY),
                     size = 0.1, colour = "grey50") +
        geom_point(size = 0.5) +
        theme_void() +
        theme(legend.position = "none")

    if (!is_discrete) {
        gg <- gg + scale_color_viridis_c()
    }

    if (label) {
        clust_data <- embedding %>%
            group_by(Cluster) %>%
            summarise(X = mean(X),
                      Y = mean(Y))

        gg <- gg +
            geom_point(data = clust_data, aes(fill = Cluster),
                       size = 10, shape = 21, colour = "white") +
            geom_text(data = clust_data, aes(label = Cluster),
                      colour = "white")
    }

    return(gg)
}


#' PAGA compare
#'
#' Compare PAGA cluster and cell embeddings
#'
#' @param clust_embedding data.frame with positions for cluster nodes
#' @param clust_edges data.frame with cluster edges
#' @param clust_thresh threshold for cluster connectivity
#' @param cell_embedding data.frame with positions for cells
#' @param cell_edges data.frame with cell edges
#' @param cell_thresh threshold for cell connectivity
#' @param colour column to use for cluster/cell colour
#' @param label whether to add cluster labels to cell plot
#'
#' @return ggplot object
plotPAGACompare <- function(clust_embedding, clust_edges, clust_thresh = 0,
                            cell_embedding, cell_edges, cell_thresh = 0,
                            colour = "Cluster", label = FALSE) {

    clusts <- plotPAGAClustGraph(clust_embedding, clust_edges, clust_thresh,
                                 colour)

    cells <- plotPAGACellGraph(cell_embedding, cell_edges, cell_thresh,
                               colour, label)

    cowplot::plot_grid(clusts, cells, nrow = 1)
}
