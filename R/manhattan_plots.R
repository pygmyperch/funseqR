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
#' @param label_type Character. Type of labels for enriched loci: "go_term", "gene_name", "uniprot_accession", or "position". Default is "go_term"
#' @param label_top_hits Integer. Number of top statistical hits to label if no functional loci. Default is 0
#' @param numeric_x_labels Logical. Use numeric labels (1,2,3,...,U) instead of chromosome names (LG1,LG2,...,U). Default is FALSE
#' @param enriched_point_size Numeric. Size of enriched loci points. Default is point_size * 1.5
#' @param enriched_point_shape Integer. Shape (pch) for enriched loci points. Default is 19 (filled circle)
#' @param enriched_point_color Character. Color for enriched loci points. Default is highlight_color
#' @param use_label_lines Logical. Use indicator lines from labels to points. Default is TRUE
#' @param signif_line_color Character. Color for significance threshold line. Default is "red"
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return ggplot2 object
#'
#' @details
#' This function creates a Manhattan plot with functional annotation highlighting:
#'
#' \\strong{Chromosome Consolidation:}
#' Uses consolidated chromosome names from the database metadata (set via define_chromosomes()).
#' Main chromosomes (e.g., LG1-LG24) are displayed separately, while scaffolds are grouped as "U".
#'
#' \\strong{Functional Highlighting:}
#' Points corresponding to functionally enriched loci (from functional_summary$loci_summary)
#' are highlighted in a different color and automatically labeled.
#'
#' \\strong{Labeling Options:}
#' Functional loci can be labeled with GO terms, gene names, UniProt accessions, or positions
#' based on the label_type parameter. Labels are automatically applied to all enriched loci.
#'
#' \\strong{Visual Features:}
#' - Alternating chromosome colors for easy visualization
#' - Automatic chromosome ordering (LG1, LG2, ..., LG24, U)
#' - Equal chromosome widths for balanced visualization
#' - Optional significance threshold line
#' - Proper x-axis spacing and labeling
#'
#' \\strong{Chromosome Width:}
#' Each chromosome gets equal width on the x-axis regardless of variant count,
#' providing balanced visualization across all chromosomes.
#'
#' @examples
#' \dontrun{
#' con <- connect_funseq_db("analysis.db")
#'
#' # Define main chromosomes for consolidation
#' define_chromosomes(con, c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8",
#'                          "LG9", "LG10", "LG11", "LG12", "LG13", "LG14", "LG15",
#'                          "LG16", "LG17", "LG18", "LG19", "LG20", "LG21", "LG22",
#'                          "LG23", "LG24"))
#'
#' # Get functional summary from enrichment analysis
#' enrich_summary <- summarize_functional_loci(con,
#'                                           enrichment_results,
#'                                           candidate_file_id,
#'                                           blast_param_id = 1)
#'
#' # Create Manhattan plot with GO term labels (default)
#' manhattan_plot <- create_functional_manhattan_plot(
#'   con,
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   functional_summary = enrich_summary,
#'   y_label = "RDA q-value",
#'   plot_title = "RDA Analysis with Functional Annotation"
#' )
#'
#' # Create Manhattan plot with gene name labels and numeric x-axis
#' manhattan_plot_genes <- create_functional_manhattan_plot(
#'   con,
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   functional_summary = enrich_summary,
#'   label_type = "gene_name",
#'   numeric_x_labels = TRUE,
#'   y_label = "RDA q-value"
#' )
#'
#' # Create Manhattan plot with custom enriched point styling
#' manhattan_plot_custom <- create_functional_manhattan_plot(
#'   con,
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   functional_summary = enrich_summary,
#'   enriched_point_size = 3,
#'   enriched_point_shape = 17,  # Triangle
#'   enriched_point_color = "red",
#'   use_label_lines = TRUE
#' )
#'
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
                                           label_type = "go_term",
                                           label_top_hits = 0,
                                           numeric_x_labels = FALSE,
                                           enriched_point_size = NULL,
                                           enriched_point_shape = 19,
                                           enriched_point_color = NULL,
                                           use_label_lines = TRUE,
                                           signif_line_color = "red",
                                           verbose = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for Manhattan plots")
  }

  if (verbose) message("Creating functional Manhattan plot...")

  # Set default values for enriched point styling
  if (is.null(enriched_point_size)) {
    enriched_point_size <- point_size * 1.5
  }
  if (is.null(enriched_point_color)) {
    enriched_point_color <- highlight_color
  }

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

  # Use consolidated chromosome names from database metadata
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

  # Handle different chromosome naming conventions
  # Extract numeric parts for proper ordering
  numeric_chrs <- unique_chrs[!unique_chrs %in% c("U", "X", "Y", "MT")]

  # For LG chromosomes, extract number after LG
  if (any(grepl("^LG", numeric_chrs))) {
    lg_numbers <- as.numeric(gsub("^LG", "", numeric_chrs[grepl("^LG", numeric_chrs)]))
    lg_numbers <- lg_numbers[!is.na(lg_numbers)]
    numeric_chrs <- paste0("LG", sort(lg_numbers))
  } else {
    # For pure numeric chromosomes
    pure_numeric <- suppressWarnings(as.numeric(numeric_chrs))
    numeric_chrs <- numeric_chrs[order(pure_numeric, na.last = TRUE)]
  }

  # Final chromosome order: numeric, then special, then unmapped
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

  # Prepare for plotting - calculate x-axis positions with equal chromosome widths
  manhattan_data <- manhattan_data[order(manhattan_data$chromosome, manhattan_data$position), ]

  # Equal width for each chromosome regardless of variant count
  unique_chrs_ordered <- levels(manhattan_data$chromosome)
  n_chrs <- length(unique_chrs_ordered)
  chr_width <- 1.0  # Each chromosome gets width of 1
  chr_spacing <- 0.2  # Gap between chromosomes

  if (verbose) message("  - Using equal chromosome widths (", n_chrs, " chromosomes)")

  # Calculate x positions with equal spacing
  manhattan_data$x_pos <- NA
  chr_midpoints <- data.frame(
    chromosome = unique_chrs_ordered,
    x_pos = NA,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(unique_chrs_ordered)) {
    chr <- unique_chrs_ordered[i]
    chr_data <- manhattan_data[manhattan_data$chromosome == chr, ]

    if (nrow(chr_data) > 0) {
      # Calculate start position for this chromosome
      chr_start <- (i - 1) * (chr_width + chr_spacing)
      chr_end <- chr_start + chr_width
      chr_midpoints$x_pos[i] <- chr_start + chr_width / 2

      # Distribute variants evenly within chromosome width
      if (nrow(chr_data) == 1) {
        manhattan_data[manhattan_data$chromosome == chr, "x_pos"] <- chr_midpoints$x_pos[i]
      } else {
        # Scale positions within chromosome width
        pos_range <- range(chr_data$position)
        if (pos_range[1] == pos_range[2]) {
          # All positions the same
          manhattan_data[manhattan_data$chromosome == chr, "x_pos"] <- chr_midpoints$x_pos[i]
        } else {
          # Scale positions proportionally within chromosome width
          scaled_pos <- chr_start + (chr_data$position - pos_range[1]) / (pos_range[2] - pos_range[1]) * chr_width
          manhattan_data[manhattan_data$chromosome == chr, "x_pos"] <- scaled_pos
        }
      }
    }
  }

  # Create numeric labels if requested
  if (numeric_x_labels) {
    # Create mapping from chromosome names to numbers
    unique_chrs_ordered <- levels(manhattan_data$chromosome)
    numeric_labels <- character(length(unique_chrs_ordered))

    numeric_counter <- 1
    for (i in seq_along(unique_chrs_ordered)) {
      chr <- unique_chrs_ordered[i]
      if (chr == "U") {
        numeric_labels[i] <- "U"
      } else {
        numeric_labels[i] <- as.character(numeric_counter)
        numeric_counter <- numeric_counter + 1
      }
    }

    # Update chr_midpoints with numeric labels
    chr_midpoints$numeric_label <- numeric_labels[match(chr_midpoints$chromosome, unique_chrs_ordered)]
    x_axis_labels <- chr_midpoints$numeric_label

    if (verbose) {
      message("  - Using numeric x-axis labels: ", paste(head(x_axis_labels, 10), collapse = ", "))
      if (length(x_axis_labels) > 10) message("    ... and ", length(x_axis_labels) - 10, " more")
    }
  } else {
    x_axis_labels <- chr_midpoints$chromosome
  }

  # Identify loci for labeling - prioritize functional loci
  manhattan_data$label <- ""

  # First, label functional loci based on label_type
  functional_loci <- which(manhattan_data$functional)
  if (length(functional_loci) > 0) {
    if (verbose) message("  - Labeling ", length(functional_loci), " functionally enriched loci")

    for (idx in functional_loci) {
      label_text <- ""

      # Find corresponding loci_summary entry for this position
      if (!is.null(functional_summary) && !is.null(functional_summary$loci_summary)) {
        loci_match <- which(functional_summary$loci_summary$chromosome == manhattan_data$chromosome[idx] &
                           functional_summary$loci_summary$position == manhattan_data$position[idx])

        if (length(loci_match) > 0) {
          loci_info <- functional_summary$loci_summary[loci_match[1], ]

          if (label_type == "go_term") {
            # Use first enriched term
            terms <- strsplit(manhattan_data$enriched_terms[idx], ";")[[1]]
            if (length(terms) > 0) {
              label_text <- trimws(terms[1])
              # Truncate if too long
              if (nchar(label_text) > 30) {
                label_text <- paste0(substr(label_text, 1, 27), "...")
              }
            }
          } else if (label_type == "gene_name" && !is.null(loci_info$gene_names) && !is.na(loci_info$gene_names) && loci_info$gene_names != "") {
            # Use gene names
            genes <- strsplit(loci_info$gene_names, ";")[[1]]
            label_text <- trimws(genes[1])
          } else if (label_type == "uniprot_accession" && !is.null(loci_info$uniprot_accession) && !is.na(loci_info$uniprot_accession)) {
            # Use UniProt accession
            label_text <- loci_info$uniprot_accession
          }
        }
      }

      # Fallback to position if no specific label found
      if (label_text == "") {
        label_text <- paste0(manhattan_data$chromosome[idx], ":", format(manhattan_data$position[idx], big.mark = ","))
      }

      manhattan_data$label[idx] <- label_text
    }
  }

  # Optionally label top statistical hits if no functional loci or if requested
  if (label_top_hits > 0 && length(functional_loci) == 0) {
    if (verbose) message("  - No functional loci found, labeling ", label_top_hits, " top statistical hits")

    if (transform_y == "neg_log10") {
      top_indices <- order(manhattan_data$y_transformed, decreasing = TRUE)[1:min(label_top_hits, nrow(manhattan_data))]
    } else {
      top_indices <- order(manhattan_data$y_value, decreasing = TRUE)[1:min(label_top_hits, nrow(manhattan_data))]
    }

    for (idx in top_indices) {
      manhattan_data$label[idx] <- paste0(manhattan_data$chromosome[idx], ":",
                                         format(manhattan_data$position[idx], big.mark = ","))
    }
  }

  # Create the plot
  p <- ggplot2::ggplot(manhattan_data, ggplot2::aes(x = x_pos, y = y_transformed)) +
    # Background points (non-functional)
    ggplot2::geom_point(data = manhattan_data[!manhattan_data$functional, ],
                       ggplot2::aes(color = chr_color),
                       size = point_size, alpha = 0.7) +
    # Highlighted functional points with custom styling
    ggplot2::geom_point(data = manhattan_data[manhattan_data$functional, ],
                       color = enriched_point_color,
                       size = enriched_point_size,
                       shape = enriched_point_shape,
                       alpha = 0.9) +
    # Manual color scale for chromosomes
    ggplot2::scale_color_identity() +
    # X-axis
    ggplot2::scale_x_continuous(
      breaks = chr_midpoints$x_pos,
      labels = x_axis_labels,
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
      # Remove all grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      # Add axis lines
      axis.line = ggplot2::element_line(color = "black", size = 0.5),
      # Add axis ticks
      axis.ticks = ggplot2::element_line(color = "black", size = 0.3),
      axis.ticks.length = ggplot2::unit(0.2, "cm"),
      # Text formatting
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "none",
      # Clean panel background
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
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
                                color = signif_line_color,
                                alpha = 0.7)
  }

  # Add labels for enriched loci
  labeled_data <- manhattan_data[manhattan_data$label != "", ]
  if (nrow(labeled_data) > 0) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      # Fallback to basic text labels if ggrepel not available
      if (use_label_lines) {
        warning("ggrepel package required for label lines. Using basic text labels without lines.")
      }
      p <- p + ggplot2::geom_text(data = labeled_data,
                                 ggplot2::aes(label = label),
                                 size = 3, vjust = -0.5, hjust = 0.5)
    } else {
      # Use ggrepel for better label positioning with indicator lines
      if (use_label_lines) {
        p <- p + ggrepel::geom_text_repel(data = labeled_data,
                                         ggplot2::aes(label = label),
                                         size = 3,
                                         max.overlaps = Inf,
                                         box.padding = 0.5,
                                         point.padding = 0.3,
                                         segment.color = "black",
                                         segment.size = 0.3,
                                         segment.alpha = 0.7,
                                         min.segment.length = 0,
                                         force = 2,
                                         force_pull = 0.5)
      } else {
        # No indicator lines
        p <- p + ggrepel::geom_text_repel(data = labeled_data,
                                         ggplot2::aes(label = label),
                                         size = 3,
                                         max.overlaps = Inf,
                                         box.padding = 0.3,
                                         point.padding = 0.3,
                                         segment.size = 0)  # No lines
      }
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
#' Creates a basic Manhattan plot without functional highlighting. Uses consolidated
#' chromosome names from database if available.
#'
#' @param con Database connection object
#' @param y_values Numeric vector of values to plot on y-axis (same order as VCF file variants)
#' @param vcf_file_id Integer. File ID of the VCF file
#' @param y_label Character. Label for y-axis. Default is "Statistical Value"
#' @param plot_title Character. Title for the plot. Default is "Manhattan Plot"
#' @param signif_threshold Numeric. Significance threshold line to draw. Default is 0.01
#' @param transform_y Character. Transform y-values: "none", "neg_log10", or "log10". Default is "neg_log10"
#' @param chr_colors Character vector. Two colors for alternating chromosomes. Default is snapper colors
#' @param point_size Numeric. Size of points. Default is 1.2
#' @param numeric_x_labels Logical. Use numeric labels (1,2,3,...,U) instead of chromosome names. Default is FALSE
#' @param signif_line_color Character. Color for significance threshold line. Default is "red"
#' @param verbose Logical. Print progress information. Default is TRUE
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Basic Manhattan plot
#' manhattan_plot <- create_manhattan_plot(
#'   con,
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   y_label = "RDA q-value",
#'   plot_title = "RDA Analysis"
#' )
#'
#' # Manhattan plot with numeric x-axis labels
#' manhattan_plot_numeric <- create_manhattan_plot(
#'   con,
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   y_label = "RDA q-value",
#'   plot_title = "RDA Analysis",
#'   numeric_x_labels = TRUE
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
                                 numeric_x_labels = FALSE,
                                 signif_line_color = "red",
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
    label_type = "position",
    label_top_hits = 0,
    numeric_x_labels = numeric_x_labels,
    signif_line_color = signif_line_color,
    verbose = verbose
  )
}

# INTERNAL HELPER FUNCTIONS

#' Get consolidated chromosome names from database metadata
#'
#' @param con Database connection object
#' @param chromosomes Character vector of chromosome names from VCF
#' @param verbose Logical. Print progress information
#' @return Character vector of consolidated chromosome names
#' @keywords internal
.get_consolidated_chromosome_names <- function(con, chromosomes, verbose = FALSE) {

  if (verbose) message("  - Using consolidated chromosome names from database metadata...")

  # Get main chromosomes from metadata
  tryCatch({
    # Check if metadata table exists
    tables <- DBI::dbListTables(con)
    if (!"metadata" %in% tables) {
      if (verbose) message("  - No metadata table found, using original names")
      return(chromosomes)
    }

    # Get main chromosomes from metadata
    result <- DBI::dbGetQuery(con, "
      SELECT value FROM metadata WHERE key = 'main_chromosomes'
    ")

    if (nrow(result) == 0) {
      if (verbose) message("  - No main chromosomes defined, using original names")
      return(chromosomes)
    }

    # Parse JSON
    main_chromosomes <- jsonlite::fromJSON(result$value[1])

    if (!is.character(main_chromosomes) || length(main_chromosomes) == 0) {
      if (verbose) message("  - Invalid main chromosomes data, using original names")
      return(chromosomes)
    }

    # Apply consolidation: main chromosomes stay as-is, others become "U"
    consolidated <- chromosomes
    consolidated[!chromosomes %in% main_chromosomes] <- "U"

    if (verbose) {
      main_count <- sum(chromosomes %in% main_chromosomes)
      scaffold_count <- sum(!chromosomes %in% main_chromosomes)
      total_count <- length(chromosomes)

      message("  - Main chromosomes (", length(main_chromosomes), "): ", paste(main_chromosomes, collapse = ", "))
      message("  - Applied consolidation to ", total_count, " variants:")
      message("    - Main chromosomes: ", main_count, " variants")
      message("    - Scaffolds -> U: ", scaffold_count, " variants")

      # Show some example mappings
      unique_original <- unique(chromosomes[!chromosomes %in% main_chromosomes])
      if (length(unique_original) > 0) {
        sample_scaffolds <- head(unique_original, 5)
        cat("  - Example scaffold mappings: ", paste(sample_scaffolds, "-> U", collapse = ", "))
        if (length(unique_original) > 5) {
          cat(" (and ", length(unique_original) - 5, " more)")
        }
        cat("\n")
      }
    }

    return(consolidated)

  }, error = function(e) {
    if (verbose) message("  - Error accessing metadata: ", e$message, ". Using original names.")
    return(chromosomes)
  })
}
