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
