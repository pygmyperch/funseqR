#' Manhattan Plot Functions for funseqR
#'
#' Functions to create Manhattan plots for genome-wide association and environmental association studies,
#' with functional annotation highlighting from GO enrichment results.

#' Create functional Manhattan plot with enrichment annotation
#'
#' Creates a Manhattan plot showing statistical values (p-values, q-values, etc.) across the genome
#' with highlighted points corresponding to functionally enriched loci from GO analysis.
#'
#' @param con Database connection object
#' @param y_values Numeric vector of values to plot on y-axis (same order as VCF file variants)
#' @param vcf_file_id Integer. File ID of the VCF file used in the analysis
#' @param functional_summary List. Output from summarize_functional_loci() containing loci_summary
#' @param y_label Character. Label for y-axis. Default is "Statistical Value"
#' @param plot_title Character. Title for the plot. Default is "Functional Manhattan Plot"
#' @param signif_threshold Numeric. Significance threshold line to draw. Default is 0.01
#' @param transform_y Character. Transform y-values: "none", "neg_log10", or "log10". Default is "neg_log10"
#' @param highlight_color Character. Color for functionally enriched points. Default is "#D0DBEE"
#' @param chr_colors Character vector. Two colors for alternating chromosomes. Default is snapper colors
#' @param point_size Numeric. Size of points. Default is 1.2
#' @param label_top_hits Integer. Number of top hits to label. Default is 10
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return ggplot2 object
#'
#' @details
#' This function creates a Manhattan plot in the style of the provided example, with:
#' - Alternating chromosome colors for easy visualization
#' - Highlighted points for functionally enriched loci
#' - Optional significance threshold line
#' - Labels for top significant hits
#' - Proper handling of unmapped contigs
#'
#' The function matches the y_values vector (in VCF order) with genomic coordinates
#' from the database and highlights any positions that appear in the functional
#' summary from GO enrichment analysis.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#' 
#' # Get functional summary from enrichment analysis
#' enrich_summary <- summarize_functional_loci(con, 
#'                                           enrichment_results, 
#'                                           candidate_file_id,
#'                                           blast_param_id = 1)
#' 
#' # Create Manhattan plot with RDA q-values
#' manhattan_plot <- create_functional_manhattan_plot(
#'   con, 
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   functional_summary = enrich_summary,
#'   y_label = "RDA q-value",
#'   plot_title = "RDA Analysis with Functional Annotation",
#'   signif_threshold = 0.01
#' )
#' print(manhattan_plot)
#' }
#'
#' @export
create_functional_manhattan_plot <- function(con, y_values, vcf_file_id, functional_summary,
                                           y_label = "Statistical Value", 
                                           plot_title = "Functional Manhattan Plot",
                                           signif_threshold = 0.01,
                                           transform_y = "neg_log10",
                                           highlight_color = "#D0DBEE",
                                           chr_colors = c("#A1B1CC", "#0E7EC0"),
                                           point_size = 1.2,
                                           label_top_hits = 10,
                                           verbose = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for Manhattan plots")
  }
  
  if (verbose) message("Creating functional Manhattan plot...")
  
  # Get genomic coordinates from database in VCF order
  vcf_coords <- DBI::dbGetQuery(con, "
    SELECT vcf_id, chromosome, position, ref, alt
    FROM vcf_data 
    WHERE file_id = ?
    ORDER BY vcf_id
  ", list(vcf_file_id))
  
  if (nrow(vcf_coords) == 0) {
    stop("No VCF data found for file_id: ", vcf_file_id)
  }
  
  # Validate y_values length
  if (length(y_values) != nrow(vcf_coords)) {
    stop("Length mismatch: y_values has ", length(y_values), 
         " values but VCF has ", nrow(vcf_coords), " variants")
  }
  
  if (verbose) message("  - Processing ", nrow(vcf_coords), " variants")
  
  # Create base Manhattan data
  manhattan_data <- data.frame(
    vcf_id = vcf_coords$vcf_id,
    chromosome = vcf_coords$chromosome,
    position = vcf_coords$position,
    ref = vcf_coords$ref,
    alt = vcf_coords$alt,
    y_value = y_values,
    stringsAsFactors = FALSE
  )
  
  # Use consolidated chromosome names from database if available
  manhattan_data$chromosome <- .get_consolidated_chromosome_names(con, manhattan_data$chromosome, verbose = verbose)
  
  # Transform y-values if requested
  if (transform_y == "neg_log10") {
    manhattan_data$y_transformed <- -log10(pmax(manhattan_data$y_value, 1e-300))  # Avoid log(0)
    if (is.null(y_label) || y_label == "Statistical Value") {
      y_label <- paste0("-log10(", y_label, ")")
    }
  } else if (transform_y == "log10") {
    manhattan_data$y_transformed <- log10(pmax(manhattan_data$y_value, 1e-300))
    if (is.null(y_label) || y_label == "Statistical Value") {
      y_label <- paste0("log10(", y_label, ")")
    }
  } else {
    manhattan_data$y_transformed <- manhattan_data$y_value
  }
  
  # Create chromosome factor with proper ordering
  unique_chrs <- unique(manhattan_data$chromosome)
  numeric_chrs <- unique_chrs[!unique_chrs %in% c("U", "X", "Y", "MT")]
  numeric_chrs <- numeric_chrs[order(as.numeric(numeric_chrs))]
  chr_levels <- c(numeric_chrs, "X", "Y", "MT", "U")
  chr_levels <- chr_levels[chr_levels %in% unique_chrs]
  
  manhattan_data$chromosome <- factor(manhattan_data$chromosome, levels = chr_levels)
  
  # Add alternating colors for chromosomes
  manhattan_data$chr_color <- ifelse(as.numeric(manhattan_data$chromosome) %% 2 == 1, 
                                   chr_colors[1], chr_colors[2])
  
  # Mark functionally enriched loci
  manhattan_data$functional <- FALSE
  manhattan_data$enriched_terms <- ""
  
  if (!is.null(functional_summary) && !is.null(functional_summary$loci_summary)) {
    loci_summary <- functional_summary$loci_summary
    
    if (nrow(loci_summary) > 0) {
      if (verbose) message("  - Highlighting ", nrow(loci_summary), " functionally enriched loci")
      
      # Match functional loci to Manhattan data
      for (i in 1:nrow(loci_summary)) {
        matches <- which(manhattan_data$chromosome == loci_summary$chromosome[i] & 
                        manhattan_data$position == loci_summary$position[i])
        
        if (length(matches) > 0) {
          manhattan_data$functional[matches] <- TRUE
          # Truncate long term lists for labeling
          terms <- loci_summary$enriched_terms[i]
          if (nchar(terms) > 100) {
            terms <- paste0(substr(terms, 1, 97), "...")
          }
          manhattan_data$enriched_terms[matches] <- terms
        }
      }
      
      func_count <- sum(manhattan_data$functional)
      if (verbose) message("  - Successfully matched ", func_count, " enriched loci")
    }
  }
  
  # Prepare for plotting - calculate cumulative positions for x-axis
  manhattan_data <- manhattan_data[order(manhattan_data$chromosome, manhattan_data$position), ]
  
  # Calculate cumulative positions for continuous x-axis
  chr_lengths <- aggregate(position ~ chromosome, data = manhattan_data, FUN = max)
  chr_lengths$cum_length <- cumsum(c(0, chr_lengths$position[-nrow(chr_lengths)]))
  
  # Add cumulative positions to data
  manhattan_data <- merge(manhattan_data, chr_lengths[, c("chromosome", "cum_length")], 
                         by = "chromosome", all.x = TRUE)
  manhattan_data$x_pos <- manhattan_data$position + manhattan_data$cum_length
  
  # Calculate chromosome midpoints for x-axis labels
  chr_midpoints <- aggregate(x_pos ~ chromosome, data = manhattan_data, FUN = function(x) mean(range(x)))
  
  # Identify top hits for labeling
  manhattan_data$label <- ""
  if (label_top_hits > 0) {
    if (transform_y == "neg_log10") {
      top_indices <- order(manhattan_data$y_transformed, decreasing = TRUE)[1:min(label_top_hits, nrow(manhattan_data))]
    } else {
      top_indices <- order(manhattan_data$y_value, decreasing = TRUE)[1:min(label_top_hits, nrow(manhattan_data))]
    }
    
    for (idx in top_indices) {
      if (manhattan_data$functional[idx]) {
        # For functional loci, use first GO term
        terms <- strsplit(manhattan_data$enriched_terms[idx], ";")[[1]]
        manhattan_data$label[idx] <- trimws(terms[1])
      } else {
        # For non-functional loci, use position
        manhattan_data$label[idx] <- paste0(manhattan_data$chromosome[idx], ":", 
                                          format(manhattan_data$position[idx], big.mark = ","))
      }
    }
  }
  
  # Create the plot
  p <- ggplot2::ggplot(manhattan_data, ggplot2::aes(x = x_pos, y = y_transformed)) +
    # Background points (non-functional)
    ggplot2::geom_point(data = manhattan_data[!manhattan_data$functional, ],
                       ggplot2::aes(color = chr_color), 
                       size = point_size, alpha = 0.7) +
    # Highlighted functional points
    ggplot2::geom_point(data = manhattan_data[manhattan_data$functional, ],
                       color = highlight_color, size = point_size * 1.5, alpha = 0.9) +
    # Manual color scale for chromosomes
    ggplot2::scale_color_identity() +
    # X-axis
    ggplot2::scale_x_continuous(
      breaks = chr_midpoints$x_pos,
      labels = chr_midpoints$chromosome,
      expand = c(0.01, 0)
    ) +
    # Y-axis
    ggplot2::scale_y_continuous(
      expand = c(0.02, 0)
    ) +
    # Labels
    ggplot2::labs(
      title = plot_title,
      x = "Chromosome",
      y = y_label
    ) +
    # Theme
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "none"
    )
  
  # Add significance threshold line
  if (!is.null(signif_threshold) && signif_threshold > 0) {
    if (transform_y == "neg_log10") {
      threshold_y <- -log10(signif_threshold)
    } else if (transform_y == "log10") {
      threshold_y <- log10(signif_threshold)
    } else {
      threshold_y <- signif_threshold
    }
    
    p <- p + ggplot2::geom_hline(yintercept = threshold_y, 
                                linetype = "dashed", 
                                color = "red", 
                                alpha = 0.7)
  }
  
  # Add labels for top hits
  labeled_data <- manhattan_data[manhattan_data$label != "", ]
  if (nrow(labeled_data) > 0) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      # Fallback to basic text labels if ggrepel not available
      p <- p + ggplot2::geom_text(data = labeled_data,
                                 ggplot2::aes(label = label),
                                 size = 3, vjust = -0.5, hjust = 0.5)
    } else {
      # Use ggrepel for better label positioning
      p <- p + ggrepel::geom_text_repel(data = labeled_data,
                                       ggplot2::aes(label = label),
                                       size = 3, 
                                       max.overlaps = Inf,
                                       box.padding = 0.3,
                                       point.padding = 0.3)
    }
  }
  
  if (verbose) {
    message("Manhattan plot created successfully")
    if (sum(manhattan_data$functional) > 0) {
      message("  - ", sum(manhattan_data$functional), " functionally enriched loci highlighted")
    }
    if (label_top_hits > 0) {
      message("  - ", nrow(labeled_data), " top hits labeled")
    }
  }
  
  return(p)
}

#' Simple Manhattan plot without functional annotation
#'
#' Creates a basic Manhattan plot without functional highlighting.
#'
#' @param con Database connection object
#' @param y_values Numeric vector of values to plot on y-axis (same order as VCF file variants)
#' @param vcf_file_id Integer. File ID of the VCF file
#' @param y_label Character. Label for y-axis. Default is "Statistical Value"
#' @param plot_title Character. Title for the plot. Default is "Manhattan Plot"
#' @param signif_threshold Numeric. Significance threshold line to draw. Default is 0.01
#' @param transform_y Character. Transform y-values: "none", "neg_log10", or "log10". Default is "neg_log10"
#' @param chr_colors Character vector. Two colors for alternating chromosomes
#' @param point_size Numeric. Size of points. Default is 1.2
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' manhattan_plot <- create_manhattan_plot(
#'   con, 
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   y_label = "RDA q-value",
#'   plot_title = "RDA Analysis"
#' )
#' }
#'
#' @export
create_manhattan_plot <- function(con, y_values, vcf_file_id,
                                 y_label = "Statistical Value", 
                                 plot_title = "Manhattan Plot",
                                 signif_threshold = 0.01,
                                 transform_y = "neg_log10",
                                 chr_colors = c("#A1B1CC", "#0E7EC0"),
                                 point_size = 1.2,
                                 verbose = TRUE) {
  
  # Call the functional version with NULL functional summary
  create_functional_manhattan_plot(
    con = con,
    y_values = y_values,
    vcf_file_id = vcf_file_id,
    functional_summary = NULL,
    y_label = y_label,
    plot_title = plot_title,
    signif_threshold = signif_threshold,
    transform_y = transform_y,
    highlight_color = NULL,
    chr_colors = chr_colors,
    point_size = point_size,
    label_top_hits = 0,
    verbose = verbose
  )
}

# INTERNAL HELPER FUNCTIONS

#' Get consolidated chromosome names from database
#' 
#' @param con Database connection object
#' @param chromosomes Character vector of chromosome names from VCF
#' @param verbose Logical. Print progress information
#' @return Character vector of consolidated chromosome names
#' @keywords internal
.get_consolidated_chromosome_names <- function(con, chromosomes, verbose = FALSE) {
  
  if (verbose) message("  - Using consolidated chromosome names from database...")
  
  # Check if chromosome consolidation table exists
  tables <- DBI::dbListTables(con)
  if (!"chromosome_consolidation" %in% tables) {
    if (verbose) message("  - No chromosome consolidation found, using original names")
    return(chromosomes)
  }
  
  # Get consolidation mapping from database
  tryCatch({
    consolidation_map <- DBI::dbGetQuery(con, "
      SELECT original_name, consolidated_name 
      FROM chromosome_consolidation
    ")
    
    if (nrow(consolidation_map) == 0) {
      if (verbose) message("  - Empty chromosome consolidation table, using original names")
      return(chromosomes)
    }
    
    # Apply consolidation mapping
    consolidated <- chromosomes
    for (i in 1:nrow(consolidation_map)) {
      original <- consolidation_map$original_name[i]
      consolidated_name <- consolidation_map$consolidated_name[i]
      consolidated[chromosomes == original] <- consolidated_name
    }
    
    if (verbose) {
      mapped_count <- sum(chromosomes %in% consolidation_map$original_name)
      total_count <- length(chromosomes)
      message("  - Applied consolidation to ", mapped_count, " of ", total_count, " variants")
      
      # Show unique mappings
      unique_mappings <- unique(data.frame(
        original = chromosomes[chromosomes %in% consolidation_map$original_name],
        consolidated = consolidated[chromosomes %in% consolidation_map$original_name],
        stringsAsFactors = FALSE
      ))
      
      if (nrow(unique_mappings) <= 30) {  # Only show if not too many
        cat("  - Chromosome mappings applied:\n")
        for (i in 1:nrow(unique_mappings)) {
          cat("    ", unique_mappings$original[i], " -> ", unique_mappings$consolidated[i], "\n")
        }
      }
    }
    
    return(consolidated)
    
  }, error = function(e) {
    if (verbose) message("  - Error accessing consolidation table: ", e$message, ". Using original names.")
    return(chromosomes)
  })
}