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
#' @param y_values Numeric vector of values to plot on y-axis (same order as VCF file variants). If NULL, uses stored statistics from database
#' @param vcf_file_id Integer. File ID of the VCF file used in the analysis
#' @param enrichment_data Data.frame. Output from compile_funseq_results() with enrichment stage. If NULL, creates basic Manhattan plot. Default is NULL
#' @param y_label Character. Label for y-axis. Default is "Statistical Value"
#' @param signif_threshold Numeric. Significance threshold line to draw. Default is 0.01
#' @param transform_y Character. Transform y-values: "none", "neg_log10", or "log10". Default is "neg_log10"
#' @param highlight_color Character. Default color for enriched points (used if enriched_point_color is NULL). Default is "#D0DBEE"
#' @param chr_colors Character vector. Two colors for alternating chromosomes. Default is snapper colors
#' @param point_size Numeric. Size of points. Default is 1.2
#' @param label_type Character. Type of labels for enriched loci: "go_term", "go_id", "gene_name", "uniprot_accession", or "position". Default is "go_term"
#' @param label_cex Numeric. Size of text labels. Default is 0.8
#' @param label_top_candidates Integer. Number of top candidate loci to label based on y_values. Uses functional annotations when available, falls back to position. Default is 0
#' @param numeric_x_labels Logical. Use numeric labels (1,2,3,...,U) instead of chromosome names (LG1,LG2,...,U). Default is FALSE
#' @param enriched_point_size Numeric. Size of enriched loci points. Default is point_size * 1.5
#' @param enriched_point_shape Integer. Shape (pch) for enriched loci points. Default is 17 (triangle)
#' @param enriched_point_color Character. Color for enriched loci points. Default is "red"
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
#' Points corresponding to functionally enriched loci (from enrichment_data with enriched = TRUE)
#' are highlighted in a different color and automatically labeled. Use dataset_type to control
#' whether to highlight candidate loci, background loci, or both.
#'
#' \\strong{Labeling Options:}
#' Enriched loci can be labeled with different annotation types using the label_type parameter:
#' \\itemize{
#'   \\item \\strong{go_term}: Enriched GO term names (e.g., "P:hemopoiesis") - human readable
#'   \\item \\strong{go_id}: Enriched GO term IDs (e.g., "GO:0030097") - compact format
#'   \\item \\strong{gene_name}: Gene names (e.g., "zfpm1") - protein/gene identifiers (displayed in italics)
#'   \\item \\strong{uniprot_accession}: UniProt accessions - database identifiers
#'   \\item \\strong{position}: Genomic position (e.g., "LG4:3814415") - coordinate fallback
#' }
#' Labels use only significantly enriched terms, not all functional annotations.
#' Gene names are automatically formatted in italics following scientific convention.
#'
#' \\strong{Top Candidate Labeling:}
#' The \\code{label_top_candidates} parameter labels the highest y_value loci using
#' intelligent per-locus fallback logic:
#' \\enumerate{
#'   \\item Try functional annotation (based on label_type: go_term, gene_name, etc.)
#'   \\item Fall back to position labels (e.g., "LG4:3,814,415") if no annotation found
#' }
#' This creates mixed labeling where some top candidates show functional information
#' while others show positions, depending on annotation availability.
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
#' # Complete workflow: annotations -> enrichment -> manhattan plot
#'
#' # Step 1: Generate annotation data
#' results <- compile_funseq_results(
#'   con,
#'   stage = "annotations",
#'   include = c("GO", "KEGG"),
#'   candidate_loci = "candidates.vcf"
#' )
#'
#' # Step 2: Run ORA analysis
#' ORA_results <- run_ORA(con, "candidates.vcf", annotation_type = "GO")
#'
#' # Step 3: Add enrichment results
#' enrichment_results <- compile_funseq_results(
#'   con,
#'   stage = "enrichment",
#'   data = results,
#'   analysis_ids = c(1, 2, 3)
#' )
#'
#' # Create Manhattan plot with enriched loci highlighted (candidate loci only)
#' manhattan_plot <- create_functional_manhattan_plot(
#'   con,
#'   y_values = my_statistical_values,  # Your p-values, q-values, etc.
#'   vcf_file_id = 1,
#'   enrichment_data = enrichment_results,
#'   y_label = "Statistical Value"
#' )
#' 
#' # Use stored statistics (after define_locus_statistics())
#' manhattan_plot_stored <- create_functional_manhattan_plot(
#'   con,
#'   vcf_file_id = 1,  # y_values = NULL uses stored statistics
#'   enrichment_data = enrichment_results,
#'   y_label = "P-value"
#' )
#'
#' # Different labeling options for enriched loci
#' # Use GO term names (default - descriptive but longer)
#' plot_go_names <- create_functional_manhattan_plot(
#'   con, y_values = my_statistical_values, vcf_file_id = 1,
#'   enrichment_data = enrichment_results, label_type = "go_term"
#' )
#'
#' # Use GO term IDs (compact format for cleaner plots)
#' plot_go_ids <- create_functional_manhattan_plot(
#'   con, y_values = my_statistical_values, vcf_file_id = 1,
#'   enrichment_data = enrichment_results, label_type = "go_id"
#' )
#'
#' # Use gene names (intermediate length)
#' plot_genes <- create_functional_manhattan_plot(
#'   con, y_values = my_statistical_values, vcf_file_id = 1,
#'   enrichment_data = enrichment_results, label_type = "gene_name"
#' )
#'
#' # Label top 10 candidates with mixed functional/position labels
#' plot_top_candidates <- create_functional_manhattan_plot(
#'   con, y_values = my_statistical_values, vcf_file_id = 1,
#'   enrichment_data = enrichment_results,
#'   label_top_candidates = 10, label_type = "gene_name"
#'   # Results in mixed labels: "zfpm1", "gata1", "LG10:28,936,085", "tbx5", etc.
#' )
#'
#' print(manhattan_plot)
#' }
#'
#' @export
create_functional_manhattan_plot <- function(con, y_values = NULL, vcf_file_id, enrichment_data = NULL,
                                           y_label = "Statistical Value",
                                           signif_threshold = 0.01,
                                           transform_y = "neg_log10",
                                           highlight_color = "#D0DBEE",
                                           chr_colors = c("#A1B1CC", "#0E7EC0"),
                                           point_size = 1.2,
                                           label_type = "go_term",
                                           label_cex = 0.8,
                                           label_top_candidates = 0,
                                           numeric_x_labels = FALSE,
                                           enriched_point_size = NULL,
                                           enriched_point_shape = 17,
                                           enriched_point_color = "red",
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

  # Handle y_values - either provided or retrieve from stored statistics
  if (is.null(y_values)) {
    if (verbose) message("  - Retrieving stored statistics from database...")
    
    # Check if stored statistics exist
    stored_stats <- DBI::dbGetQuery(con, "SELECT COUNT(*) as count FROM locus_statistics")$count
    if (stored_stats == 0) {
      stop("No y_values provided and no stored statistics found. Use define_locus_statistics() first or provide y_values.")
    }
    
    # Retrieve statistics matching VCF coordinates
    stats_query <- DBI::dbGetQuery(con, "
      SELECT v.vcf_id, v.chromosome, v.position, COALESCE(s.statistic, NA) as statistic
      FROM vcf_data v
      LEFT JOIN locus_statistics s ON (v.chromosome = s.chromosome AND v.position = s.position)
      WHERE v.file_id = ?
      ORDER BY v.vcf_id
    ", list(vcf_file_id))
    
    # Check for missing statistics
    missing_stats <- sum(is.na(stats_query$statistic))
    if (missing_stats > 0) {
      if (missing_stats == nrow(stats_query)) {
        stop("No stored statistics match VCF coordinates. Ensure define_locus_statistics() used same coordinate system.")
      } else {
        warning("Missing statistics for ", missing_stats, " out of ", nrow(stats_query), " VCF variants. Using NA values.")
      }
    }
    
    y_values <- stats_query$statistic
    if (verbose) message("  - Retrieved ", sum(!is.na(y_values)), " stored statistics (", missing_stats, " missing)")
    
  } else {
    # Validate provided y_values length
    if (length(y_values) != nrow(vcf_coords)) {
      stop("Length mismatch: y_values has ", length(y_values),
           " values but VCF has ", nrow(vcf_coords), " variants")
    }
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
    y_label <- bquote(-log[10](.(y_label)))
  } else if (transform_y == "log10") {
    manhattan_data$y_transformed <- log10(pmax(manhattan_data$y_value, 1e-300))
    y_label <- bquote(log[10](.(y_label)))
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

  if (!is.null(enrichment_data)) {
    # Filter enriched candidate loci
    enriched_loci <- enrichment_data[enrichment_data$enriched == TRUE &
                                   enrichment_data$dataset_type == "candidate", ]

    if (nrow(enriched_loci) > 0) {
      if (verbose) message("  - Highlighting ", nrow(enriched_loci), " functionally enriched candidate loci")

      # Match enriched loci to Manhattan data
      for (i in 1:nrow(enriched_loci)) {
        matches <- which(manhattan_data$chromosome == enriched_loci$chromosome[i] &
                        manhattan_data$position == enriched_loci$position[i])

        if (length(matches) > 0) {
          manhattan_data$functional[matches] <- TRUE
          # Use enriched_terms for labeling (truncate if long)
          terms <- enriched_loci$enriched_terms[i]
          if (!is.na(terms) && terms != "") {
            if (nchar(terms) > 100) {
              terms <- paste0(substr(terms, 1, 97), "...")
            }
            manhattan_data$enriched_terms[matches] <- terms
          }
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
  manhattan_data$is_gene_name <- FALSE  # Track which labels are gene names for italic formatting

  # First, label functional loci based on label_type
  functional_loci <- which(manhattan_data$functional)
  if (length(functional_loci) > 0) {
    if (verbose) message("  - Labeling ", length(functional_loci), " functionally enriched loci")

    for (idx in functional_loci) {
      label_text <- ""

      # Find corresponding enrichment data entry for this position
      if (!is.null(enrichment_data)) {
        # Filter enriched candidate loci
        loci_data <- enrichment_data[enrichment_data$enriched == TRUE &
                                   enrichment_data$dataset_type == "candidate", ]

        loci_match <- which(loci_data$chromosome == manhattan_data$chromosome[idx] &
                           loci_data$position == manhattan_data$position[idx])

        if (length(loci_match) > 0) {
          loci_info <- loci_data[loci_match[1], ]

          if (label_type == "go_term") {
            # Use first enriched GO term name
            if (!is.na(loci_info$enriched_terms) && loci_info$enriched_terms != "") {
              terms <- strsplit(loci_info$enriched_terms, ";")[[1]]
              if (length(terms) > 0) {
                label_text <- trimws(terms[1])
                # Truncate if too long
                if (nchar(label_text) > 30) {
                  label_text <- paste0(substr(label_text, 1, 27), "...")
                }
              }
            }
          } else if (label_type == "go_id") {
            # Use first enriched GO term ID (compact format)
            if ("enriched_term_ids" %in% names(loci_info) && !is.na(loci_info$enriched_term_ids) && loci_info$enriched_term_ids != "") {
              term_ids <- strsplit(loci_info$enriched_term_ids, ";")[[1]]
              if (length(term_ids) > 0) {
                label_text <- trimws(term_ids[1])
              }
            }
          } else if (label_type == "gene_name" && !is.null(loci_info$gene_name) && !is.na(loci_info$gene_name) && loci_info$gene_name != "") {
            # Use gene names
            genes <- strsplit(loci_info$gene_name, ";")[[1]]
            label_text <- trimws(genes[1])
            manhattan_data$is_gene_name[idx] <- TRUE  # Mark as gene name for italic formatting
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

  # Optionally label top candidate loci (works independently of enrichment)
  if (label_top_candidates > 0) {
    if (verbose) message("  - Labeling ", label_top_candidates, " top candidate loci")

    if (transform_y == "neg_log10") {
      top_indices <- order(manhattan_data$y_transformed, decreasing = TRUE)[1:min(label_top_candidates, nrow(manhattan_data))]
    } else {
      top_indices <- order(manhattan_data$y_value, decreasing = TRUE)[1:min(label_top_candidates, nrow(manhattan_data))]
    }

    for (idx in top_indices) {
      # Only label if not already labeled by enrichment
      if (manhattan_data$label[idx] == "") {
        label_text <- ""

        # Try to find functional annotation for this position
        if (!is.null(enrichment_data)) {
          # Search ALL candidate loci (not just enriched=TRUE)
          all_candidates <- enrichment_data[enrichment_data$dataset_type == "candidate", ]

          loci_match <- which(all_candidates$chromosome == manhattan_data$chromosome[idx] &
                             all_candidates$position == manhattan_data$position[idx])

          if (length(loci_match) > 0) {
            loci_info <- all_candidates[loci_match[1], ]

            # Apply same labeling logic as enriched loci
            if (label_type == "go_term") {
              # Use enriched terms if available, otherwise use general GO terms
              if ("enriched_terms" %in% names(loci_info) && !is.na(loci_info$enriched_terms) && loci_info$enriched_terms != "") {
                terms <- strsplit(loci_info$enriched_terms, ";")[[1]]
                if (length(terms) > 0) {
                  label_text <- trimws(terms[1])
                  if (nchar(label_text) > 30) {
                    label_text <- paste0(substr(label_text, 1, 27), "...")
                  }
                }
              } else if (!is.na(loci_info$go_names) && loci_info$go_names != "") {
                terms <- strsplit(loci_info$go_names, ";")[[1]]
                if (length(terms) > 0) {
                  label_text <- trimws(terms[1])
                  if (nchar(label_text) > 30) {
                    label_text <- paste0(substr(label_text, 1, 27), "...")
                  }
                }
              }
            } else if (label_type == "go_id") {
              # Use enriched term IDs if available, otherwise use general GO IDs
              if ("enriched_term_ids" %in% names(loci_info) && !is.na(loci_info$enriched_term_ids) && loci_info$enriched_term_ids != "") {
                term_ids <- strsplit(loci_info$enriched_term_ids, ";")[[1]]
                if (length(term_ids) > 0) {
                  label_text <- trimws(term_ids[1])
                }
              } else if (!is.na(loci_info$go_terms) && loci_info$go_terms != "") {
                go_ids <- strsplit(loci_info$go_terms, ";")[[1]]
                if (length(go_ids) > 0) {
                  label_text <- trimws(go_ids[1])
                }
              }
            } else if (label_type == "gene_name" && !is.na(loci_info$gene_name) && loci_info$gene_name != "") {
              # Use gene names
              genes <- strsplit(loci_info$gene_name, ";")[[1]]
              label_text <- trimws(genes[1])
              manhattan_data$is_gene_name[idx] <- TRUE  # Mark as gene name for italic formatting
            } else if (label_type == "uniprot_accession" && !is.na(loci_info$uniprot_accession) && loci_info$uniprot_accession != "") {
              # Use UniProt accession
              label_text <- loci_info$uniprot_accession
            }
          }
        }

        # Fallback to position if no functional annotation found
        if (label_text == "") {
          label_text <- paste0(manhattan_data$chromosome[idx], ":",
                              format(manhattan_data$position[idx], big.mark = ","))
        }

        manhattan_data$label[idx] <- label_text
      }
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
                                 ggplot2::aes(label = label, fontface = ifelse(is_gene_name, "italic", "plain")),
                                 size = label_cex * 3, vjust = -0.5, hjust = 0.5)
    } else {
      # Use ggrepel for better label positioning with indicator lines
      if (use_label_lines) {
        p <- p + ggrepel::geom_text_repel(data = labeled_data,
                                         ggplot2::aes(label = label, fontface = ifelse(is_gene_name, "italic", "plain")),
                                         size = label_cex * 3,
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
                                         ggplot2::aes(label = label, fontface = ifelse(is_gene_name, "italic", "plain")),
                                         size = label_cex * 3,
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
    if (label_top_candidates > 0) {
      message("  - ", nrow(labeled_data), " candidate loci labeled")
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
#' @param y_values Numeric vector of values to plot on y-axis (same order as VCF file variants). If NULL, uses stored statistics from database
#' @param vcf_file_id Integer. File ID of the VCF file
#' @param y_label Character. Label for y-axis. Default is "Statistical Value"
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
#'   y_label = "RDA q-value"
#' )
#'
#' # Manhattan plot with numeric x-axis labels
#' manhattan_plot_numeric <- create_manhattan_plot(
#'   con,
#'   y_values = rda.simple.pq$q.values,
#'   vcf_file_id = 1,
#'   y_label = "RDA q-value",
#'   numeric_x_labels = TRUE
#' )
#' }
#'
#' @export
create_manhattan_plot <- function(con, y_values = NULL, vcf_file_id,
                                 y_label = "Statistical Value",
                                 signif_threshold = 0.01,
                                 transform_y = "neg_log10",
                                 chr_colors = c("#A1B1CC", "#0E7EC0"),
                                 point_size = 1.2,
                                 numeric_x_labels = FALSE,
                                 signif_line_color = "red",
                                 verbose = TRUE) {

  # Call the functional version with NULL enrichment data
  create_functional_manhattan_plot(
    con = con,
    y_values = y_values,
    vcf_file_id = vcf_file_id,
    enrichment_data = NULL,
    y_label = y_label,
    signif_threshold = signif_threshold,
    transform_y = transform_y,
    highlight_color = NULL,
    chr_colors = chr_colors,
    point_size = point_size,
    label_type = "position",
    label_top_candidates = 0,
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
