#' Functional Genomics Visualization Functions
#'
#' Functions for creating visualizations that link functional enrichment results
#' back to genomic positions, including Manhattan plots and regional views.

# INTERNAL HELPER FUNCTIONS (NOT EXPORTED)

#' Consolidate chromosome names for plotting
#' @param chromosome_data Data frame with chromosome column
#' @param consolidation_map Optional named vector for chromosome mapping
#' @return Data frame with consolidated chromosome names and numeric positions
#' @keywords internal
.consolidate_chromosomes_for_plot <- function(chromosome_data, consolidation_map = NULL) {
  if (nrow(chromosome_data) == 0) {
    return(chromosome_data)
  }
  
  # Use provided consolidation map or create default
  if (is.null(consolidation_map)) {
    # Default consolidation: try to extract numeric parts
    unique_chroms <- unique(chromosome_data$chromosome)
    
    # Extract numeric parts and sort
    numeric_chroms <- unique_chroms[grepl("^[0-9]+$", unique_chroms)]
    if (length(numeric_chroms) > 0) {
      numeric_chroms <- numeric_chroms[order(as.numeric(numeric_chroms))]
    }
    
    # Handle sex chromosomes and other standard names
    special_chroms <- c("X", "Y", "Z", "W", "MT", "M")
    present_special <- intersect(unique_chroms, special_chroms)
    
    # Handle remaining chromosomes (scaffolds, contigs, etc.)
    remaining_chroms <- setdiff(unique_chroms, c(numeric_chroms, present_special))
    
    # Create consolidated order
    ordered_chroms <- c(numeric_chroms, present_special, sort(remaining_chroms))
    consolidation_map <- setNames(paste0("Chr", seq_along(ordered_chroms)), ordered_chroms)
  }
  
  # Apply consolidation
  chromosome_data$chromosome_consolidated <- consolidation_map[chromosome_data$chromosome]
  
  # Handle unmapped chromosomes
  unmapped <- is.na(chromosome_data$chromosome_consolidated)
  if (any(unmapped)) {
    chromosome_data$chromosome_consolidated[unmapped] <- chromosome_data$chromosome[unmapped]
  }
  
  # Create numeric chromosome positions for plotting
  unique_consolidated <- unique(chromosome_data$chromosome_consolidated)
  chrom_positions <- setNames(seq_along(unique_consolidated), unique_consolidated)
  chromosome_data$chromosome_numeric <- chrom_positions[chromosome_data$chromosome_consolidated]
  
  return(chromosome_data)
}

#' Merge genomic data with GEA results
#' @param loci_data Data frame with genomic loci from traceability analysis
#' @param gea_results Data frame with GEA results (chromosome, position, metric)
#' @param position_tolerance Integer. Tolerance for position matching (bp)
#' @return Merged data frame
#' @keywords internal
.merge_loci_with_gea <- function(loci_data, gea_results, position_tolerance = 0) {
  if (nrow(loci_data) == 0 || nrow(gea_results) == 0) {
    return(data.frame())
  }
  
  # Standardize column names
  if (!"chromosome" %in% names(gea_results)) {
    # Try common alternative names
    chrom_cols <- c("chrom", "chr", "CHROM", "CHR", "scaffold", "contig")
    chrom_col <- intersect(chrom_cols, names(gea_results))
    if (length(chrom_col) > 0) {
      names(gea_results)[names(gea_results) == chrom_col[1]] <- "chromosome"
    } else {
      stop("Could not find chromosome column in gea_results. Expected one of: ", 
           paste(c("chromosome", chrom_cols), collapse = ", "))
    }
  }
  
  if (!"position" %in% names(gea_results)) {
    pos_cols <- c("pos", "POS", "bp", "start")
    pos_col <- intersect(pos_cols, names(gea_results))
    if (length(pos_col) > 0) {
      names(gea_results)[names(gea_results) == pos_col[1]] <- "position"
    } else {
      stop("Could not find position column in gea_results. Expected one of: ", 
           paste(c("position", pos_cols), collapse = ", "))
    }
  }
  
  # Convert positions to numeric
  gea_results$position <- as.numeric(gea_results$position)
  loci_data$position <- as.numeric(loci_data$position)
  
  if (position_tolerance == 0) {
    # Exact matching
    merged <- merge(loci_data, gea_results, by = c("chromosome", "position"), all.x = TRUE)
  } else {
    # Fuzzy matching within tolerance
    merged_list <- list()
    
    for (i in seq_nrow(loci_data)) {
      locus <- loci_data[i, ]
      
      # Find GEA results within tolerance
      matching_gea <- gea_results[
        gea_results$chromosome == locus$chromosome &
        abs(gea_results$position - locus$position) <= position_tolerance,
      ]
      
      if (nrow(matching_gea) > 0) {
        # Use closest match if multiple found
        if (nrow(matching_gea) > 1) {
          distances <- abs(matching_gea$position - locus$position)
          matching_gea <- matching_gea[which.min(distances), ]
        }
        
        merged_row <- cbind(locus, matching_gea[, !names(matching_gea) %in% c("chromosome", "position")])
      } else {
        # No match found
        gea_cols <- setdiff(names(gea_results), c("chromosome", "position"))
        na_values <- setNames(rep(NA, length(gea_cols)), gea_cols)
        merged_row <- cbind(locus, t(na_values))
      }
      
      merged_list[[i]] <- merged_row
    }
    
    merged <- do.call(rbind, merged_list)
  }
  
  return(merged)
}

#' Optimize labels for Manhattan plot to avoid overcrowding
#' @param plot_data Data frame with plotting data
#' @param label_column Character. Column name containing labels
#' @param max_labels Integer. Maximum number of labels to show
#' @param priority_metric Character. Column name for prioritizing labels
#' @return Logical vector indicating which points to label
#' @keywords internal
.optimize_manhattan_labels <- function(plot_data, label_column, max_labels = 20, 
                                      priority_metric = "p_adjusted") {
  if (nrow(plot_data) == 0 || max_labels <= 0) {
    return(rep(FALSE, nrow(plot_data)))
  }
  
  # Initialize label vector
  show_label <- rep(FALSE, nrow(plot_data))
  
  # Remove points with missing labels or priority metric
  valid_points <- !is.na(plot_data[[label_column]]) & !is.na(plot_data[[priority_metric]])
  
  if (sum(valid_points) == 0) {
    return(show_label)
  }
  
  # Sort by priority metric (assume lower is better for p-values)
  valid_data <- plot_data[valid_points, ]
  if (priority_metric %in% c("p_value", "p_adjusted", "fdr")) {
    # Lower p-values have higher priority
    priority_order <- order(valid_data[[priority_metric]])
  } else {
    # Higher values have higher priority (e.g., fold enrichment)
    priority_order <- order(valid_data[[priority_metric]], decreasing = TRUE)
  }
  
  # Select top priority points up to max_labels
  n_labels <- min(max_labels, length(priority_order))
  selected_indices <- which(valid_points)[priority_order[1:n_labels]]
  show_label[selected_indices] <- TRUE
  
  return(show_label)
}

#' Create chromosome breaks and labels for Manhattan plot
#' @param plot_data Data frame with chromosome_numeric and chromosome_consolidated columns
#' @return List with breaks and labels for x-axis
#' @keywords internal
.create_chromosome_breaks <- function(plot_data) {
  if (nrow(plot_data) == 0) {
    return(list(breaks = numeric(0), labels = character(0)))
  }
  
  # Calculate midpoints for each chromosome
  chrom_summary <- aggregate(
    chromosome_numeric ~ chromosome_consolidated,
    data = plot_data,
    FUN = function(x) mean(range(x))
  )
  
  breaks <- chrom_summary$chromosome_numeric
  labels <- chrom_summary$chromosome_consolidated
  
  # Sort by numeric position
  order_idx <- order(breaks)
  
  return(list(
    breaks = breaks[order_idx],
    labels = labels[order_idx]
  ))
}

#' Apply y-axis transformation
#' @param values Numeric vector of values to transform
#' @param transform Character. Transformation type
#' @return Transformed numeric vector
#' @keywords internal
.apply_y_transform <- function(values, transform = "-log10") {
  switch(transform,
    "-log10" = -log10(pmax(values, .Machine$double.eps)),
    "log10" = log10(pmax(values, .Machine$double.eps)),
    "log2" = log2(pmax(values, .Machine$double.eps)),
    "-log2" = -log2(pmax(values, .Machine$double.eps)),
    "identity" = values,
    "sqrt" = sqrt(pmax(values, 0)),
    {
      warning("Unknown transform '", transform, "', using identity")
      values
    }
  )
}

# EXPORTED FUNCTIONS

#' Create functional Manhattan plot
#'
#' Create a Manhattan plot that displays genomic positions colored and labeled by
#' functional enrichment results, with user-provided GEA statistics on the y-axis.
#' This visualization links functional predictions back to genomic evidence.
#'
#' @param con Database connection object
#' @param enrichment_results Data frame or list of GO enrichment results
#' @param candidate_file_id Integer. File ID of candidate dataset
#' @param gea_results Data frame with GEA results containing chromosome, position, and metric columns
#' @param y_metric Character. Column name in gea_results to use for y-axis. Default is "p_value"
#' @param y_transform Character. Transformation for y-axis: "-log10", "log10", "identity", etc. Default is "-log10"
#' @param significance_threshold Numeric. FDR threshold for enrichment significance. Default is 0.05
#' @param label_terms Character. Which terms to label: "significant", "top", or "none". Default is "significant"
#' @param max_labels Integer. Maximum number of labels to show. Default is 20
#' @param consolidate_scaffolds Logical. Use consolidated chromosome names. Default is TRUE
#' @param color_by Character. Color points by: "go_category", "fold_enrichment", or "significance". Default is "go_category"
#' @param point_size_by Character. Size points by: "fold_enrichment", "significance", or NULL. Default is "fold_enrichment"
#' @param position_tolerance Integer. Tolerance for matching genomic positions (bp). Default is 0
#' @param title Character. Plot title. Default is auto-generated
#' @param theme_function Function. ggplot2 theme function. Default is theme_minimal
#'
#' @return ggplot object
#'
#' @details
#' This function creates a Manhattan plot that combines genomic positions (x-axis)
#' with user-provided GEA statistics (y-axis), while highlighting and labeling
#' points based on functional enrichment results. This enables visual validation
#' that functionally enriched loci are also those with strong environmental
#' associations.
#'
#' The function automatically handles:
#' \itemize{
#'   \item Chromosome consolidation for clean visualization
#'   \item Intelligent label placement to avoid overcrowding
#'   \item Color coding by functional categories or enrichment strength
#'   \item Point sizing by enrichment metrics
#'   \item Flexible y-axis transformations
#' }
#'
#' @examples
#' \dontrun{
#' # Basic functional Manhattan plot
#' con <- connect_funseq_db("analysis.db")
#' enrich_res <- run_go_enrichment_workflow(con, "candidates.vcf")
#' 
#' # Prepare GEA results data frame
#' gea_data <- data.frame(
#'   chromosome = c("1", "1", "2", "3"),
#'   position = c(1000000, 2000000, 1500000, 800000),
#'   p_value = c(0.001, 0.0001, 0.005, 0.05)
#' )
#' 
#' # Create Manhattan plot
#' manhattan_plot <- create_functional_manhattan_plot(
#'   con, enrich_res$enrichment_results, enrich_res$candidate_import$file_id,
#'   gea_data, y_metric = "p_value"
#' )
#' 
#' print(manhattan_plot)
#' 
#' # Customize visualization
#' custom_plot <- create_functional_manhattan_plot(
#'   con, enrich_res$enrichment_results, enrich_res$candidate_import$file_id,
#'   gea_data,
#'   y_metric = "p_value",
#'   color_by = "fold_enrichment",
#'   point_size_by = "significance",
#'   max_labels = 10,
#'   title = "Environmental Association vs Functional Enrichment"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_color_viridis_c
#' @importFrom ggplot2 scale_size_continuous labs theme_minimal theme element_text
#' @importFrom ggplot2 scale_x_continuous guides guide_legend guide_colorbar
#' @export
create_functional_manhattan_plot <- function(con, enrichment_results, candidate_file_id, gea_results,
                                            y_metric = "p_value", y_transform = "-log10",
                                            significance_threshold = 0.05,
                                            label_terms = "significant", max_labels = 20,
                                            consolidate_scaffolds = TRUE,
                                            color_by = "go_category", point_size_by = "fold_enrichment",
                                            position_tolerance = 0, title = NULL,
                                            theme_function = ggplot2::theme_minimal) {
  
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    message("ggrepel package recommended for better label placement. Install with: install.packages('ggrepel')")
  }
  
  # Validate inputs
  if (!is_db_connected(con)) {
    stop("Database connection is not valid")
  }
  
  if (!is.data.frame(gea_results)) {
    stop("gea_results must be a data frame")
  }
  
  if (!y_metric %in% names(gea_results)) {
    stop("y_metric '", y_metric, "' not found in gea_results columns: ", 
         paste(names(gea_results), collapse = ", "))
  }
  
  # Trace enriched terms to genomic loci
  loci_data <- trace_enriched_terms_to_loci(
    con, enrichment_results, candidate_file_id,
    significance_threshold = significance_threshold,
    include_blast_metrics = FALSE
  )
  
  if (nrow(loci_data) == 0) {
    stop("No enriched loci found for plotting")
  }
  
  # Merge with GEA results
  plot_data <- .merge_loci_with_gea(loci_data, gea_results, position_tolerance)
  
  # Remove points without GEA data
  plot_data <- plot_data[!is.na(plot_data[[y_metric]]), ]
  
  if (nrow(plot_data) == 0) {
    stop("No overlapping positions found between enriched loci and GEA results")
  }
  
  # Consolidate chromosomes if requested
  if (consolidate_scaffolds) {
    plot_data <- .consolidate_chromosomes_for_plot(plot_data)
  } else {
    plot_data$chromosome_consolidated <- plot_data$chromosome
    plot_data$chromosome_numeric <- as.numeric(factor(plot_data$chromosome))
  }
  
  # Apply y-axis transformation
  plot_data$y_transformed <- .apply_y_transform(plot_data[[y_metric]], y_transform)
  
  # Determine which points to label
  if (label_terms == "none") {
    plot_data$show_label <- FALSE
  } else {
    priority_metric <- if (label_terms == "significant") "p_adjusted" else "fold_enrichment"
    plot_data$show_label <- .optimize_manhattan_labels(
      plot_data, "go_term", max_labels, priority_metric
    )
  }
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = chromosome_numeric, y = y_transformed))
  
  # Add points with coloring
  if (color_by == "go_category") {
    # Color by GO category
    category_colors <- c("BP" = "#1f77b4", "MF" = "#ff7f0e", "CC" = "#2ca02c")
    
    if (is.null(point_size_by)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = go_category), alpha = 0.7, size = 2)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = go_category, size = !!rlang::sym(point_size_by)), alpha = 0.7)
    }
    
    p <- p + ggplot2::scale_color_manual(
      name = "GO Category",
      values = category_colors,
      labels = c("BP" = "Biological Process", "MF" = "Molecular Function", "CC" = "Cellular Component")
    )
    
  } else if (color_by == "fold_enrichment") {
    # Color by fold enrichment
    if (is.null(point_size_by)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = fold_enrichment), alpha = 0.7, size = 2)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = fold_enrichment, size = !!rlang::sym(point_size_by)), alpha = 0.7)
    }
    
    p <- p + ggplot2::scale_color_viridis_c(name = "Fold\nEnrichment", option = "plasma")
    
  } else if (color_by == "significance") {
    # Color by significance level
    if (is.null(point_size_by)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = -log10(p_adjusted)), alpha = 0.7, size = 2)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = -log10(p_adjusted), size = !!rlang::sym(point_size_by)), alpha = 0.7)
    }
    
    p <- p + ggplot2::scale_color_viridis_c(name = "-log10\n(FDR)", option = "viridis")
  }
  
  # Add point sizing if specified
  if (!is.null(point_size_by)) {
    size_name <- switch(point_size_by,
      "fold_enrichment" = "Fold\nEnrichment",
      "significance" = "-log10\n(FDR)",
      point_size_by
    )
    
    if (point_size_by == "significance") {
      plot_data$size_metric <- -log10(plot_data$p_adjusted)
    } else {
      plot_data$size_metric <- plot_data[[point_size_by]]
    }
    
    p <- p + ggplot2::scale_size_continuous(
      name = size_name,
      range = c(1, 5),
      guide = ggplot2::guide_legend(override.aes = list(alpha = 1))
    )
  }
  
  # Add labels if requested
  if (any(plot_data$show_label)) {
    label_data <- plot_data[plot_data$show_label, ]
    
    # Truncate long GO terms for labeling
    label_data$go_term_short <- stringr::str_trunc(label_data$go_term, 40)
    
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = go_term_short),
        size = 3, alpha = 0.8, max.overlaps = Inf,
        box.padding = 0.3, point.padding = 0.3
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(label = go_term_short),
        size = 3, alpha = 0.8, vjust = -0.5
      )
    }
  }
  
  # Set up x-axis
  chrom_breaks <- .create_chromosome_breaks(plot_data)
  p <- p + ggplot2::scale_x_continuous(
    name = "Chromosome",
    breaks = chrom_breaks$breaks,
    labels = chrom_breaks$labels
  )
  
  # Set up y-axis
  y_label <- switch(y_transform,
    "-log10" = paste0("-log10(", y_metric, ")"),
    "log10" = paste0("log10(", y_metric, ")"),
    "identity" = y_metric,
    paste0(y_transform, "(", y_metric, ")")
  )
  
  p <- p + ggplot2::labs(y = y_label)
  
  # Add title
  if (is.null(title)) {
    title <- "Functional Manhattan Plot: Environmental Association vs GO Enrichment"
  }
  p <- p + ggplot2::labs(title = title)
  
  # Apply theme
  p <- p + theme_function() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

#' Create enrichment region plot
#'
#' Create a detailed regional view of functional enrichment around specific
#' genomic loci, showing gene structure and annotation quality.
#'
#' @param con Database connection object
#' @param enrichment_results Data frame or list of GO enrichment results
#' @param candidate_file_id Integer. File ID of candidate dataset
#' @param chromosome Character. Target chromosome
#' @param start Integer. Start position for region
#' @param end Integer. End position for region
#' @param include_blast_metrics Logical. Include BLAST quality information. Default is TRUE
#'
#' @return ggplot object showing regional enrichment details
#'
#' @examples
#' \dontrun{
#' # Create regional plot for specific genomic region
#' region_plot <- create_enrichment_region_plot(
#'   con, enrichment_results, candidate_file_id,
#'   chromosome = "1", start = 1000000, end = 2000000
#' )
#' print(region_plot)
#' }
#'
#' @export
create_enrichment_region_plot <- function(con, enrichment_results, candidate_file_id,
                                         chromosome, start, end, include_blast_metrics = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  
  # Get all loci data for the region
  loci_data <- trace_enriched_terms_to_loci(
    con, enrichment_results, candidate_file_id,
    include_blast_metrics = include_blast_metrics
  )
  
  # Filter to specified region
  region_data <- loci_data[
    loci_data$chromosome == chromosome &
    loci_data$position >= start &
    loci_data$position <= end,
  ]
  
  if (nrow(region_data) == 0) {
    stop("No enriched loci found in specified region: ", chromosome, ":", start, "-", end)
  }
  
  # Create position-based plot
  p <- ggplot2::ggplot(region_data, ggplot2::aes(x = position, y = -log10(p_adjusted)))
  
  # Add points colored by GO category
  p <- p + ggplot2::geom_point(
    ggplot2::aes(color = go_category, size = fold_enrichment),
    alpha = 0.7
  )
  
  # Color by GO category
  category_colors <- c("BP" = "#1f77b4", "MF" = "#ff7f0e", "CC" = "#2ca02c")
  p <- p + ggplot2::scale_color_manual(
    name = "GO Category",
    values = category_colors,
    labels = c("BP" = "Biological Process", "MF" = "Molecular Function", "CC" = "Cellular Component")
  )
  
  # Size by fold enrichment
  p <- p + ggplot2::scale_size_continuous(
    name = "Fold\nEnrichment",
    range = c(2, 8)
  )
  
  # Add GO term labels
  p <- p + ggplot2::geom_text(
    ggplot2::aes(label = stringr::str_trunc(go_term, 30)),
    vjust = -0.5, size = 3, alpha = 0.8
  )
  
  # Formatting
  p <- p + ggplot2::labs(
    title = paste0("Functional Enrichment in Region ", chromosome, ":", 
                  format(start, big.mark = ","), "-", format(end, big.mark = ",")),
    x = paste0("Position on Chromosome ", chromosome),
    y = "-log10(FDR)"
  )
  
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}