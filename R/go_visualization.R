#' GO Enrichment Visualization Functions
#'
#' This file contains functions for creating various visualizations of GO enrichment results,
#' including bubble plots, treemaps, and interactive network diagrams.
#'

#' Create GO enrichment bubble plot
#'
#' @param enrichment_results Data frame. Output from perform_go_enrichment()
#' @param max_terms Integer. Maximum number of terms to show. Default is 20
#' @param min_fold_enrichment Numeric. Minimum fold enrichment to display. Default is 1.5
#' @param title Character. Plot title. Default is auto-generated
#' @param color_palette Character. Color palette for significance. Default is "plasma"
#' @param significance_threshold Numeric. FDR threshold used for significance. Default is 0.05
#'
#' @return ggplot object
#'
#' @details
#' Creates a bubble plot showing GO terms ordered by fold enrichment, with
#' bubble size representing gene count and color representing significance level.
#'
#' @examples
#' \dontrun{
#' bp_results <- perform_go_enrichment(go_data, "BP")
#' p1 <- create_go_bubble_plot(bp_results)
#' print(p1)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c scale_size_continuous
#' @importFrom ggplot2 scale_x_log10 labs theme_minimal theme element_text
#' @importFrom dplyr filter arrange slice_head mutate
#' @importFrom stringr str_trunc
#' @export
create_go_bubble_plot <- function(enrichment_results, max_terms = 20, min_fold_enrichment = 1.5, 
                                  title = NULL, color_palette = "plasma", significance_threshold = 0.05) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("stringr package is required for this function")
  }
  
  if (nrow(enrichment_results) == 0) {
    return(ggplot2::ggplot() + 
           ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No enriched terms to display", size = 5) +
           ggplot2::theme_void())
  }
  
  # Prepare data for plotting
  sig_results <- enrichment_results[enrichment_results$significance_level %in% c("significant", "highly_significant"), ]
  sig_results <- sig_results[sig_results$fold_enrichment >= min_fold_enrichment, ]
  sig_results <- sig_results[order(sig_results$p_adjusted), ]
  
  if (nrow(sig_results) > max_terms) {
    sig_results <- sig_results[1:max_terms, ]
  }
  
  plot_data <- sig_results
  plot_data$go_term_short <- stringr::str_trunc(plot_data$go_term, 50)
  plot_data$neg_log10_padj <- -log10(plot_data$p_adjusted)
  # Ensure minimum size for very small p-values
  plot_data$neg_log10_padj <- pmax(plot_data$neg_log10_padj, 1.3)  # -log10(0.05) = 1.3
  
  if (nrow(plot_data) == 0) {
    return(ggplot2::ggplot() + 
           ggplot2::annotate("text", x = 0.5, y = 0.5, 
                           label = paste0("No terms with fold enrichment >= ", min_fold_enrichment, " and FDR < 0.05"), 
                           size = 4) +
           ggplot2::theme_void())
  }
  
  # Auto-generate title if not provided
  if (is.null(title)) {
    ontology_names <- c("BP" = "Biological Process", "MF" = "Molecular Function", "CC" = "Cellular Component")
    ontology_name <- ontology_names[plot_data$go_category[1]]
    if (is.na(ontology_name)) ontology_name <- plot_data$go_category[1]
    title <- paste0("GO Enrichment: ", ontology_name)
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fold_enrichment, y = reorder(go_term_short, fold_enrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = foreground_count, color = neg_log10_padj), alpha = 0.7) +
    ggplot2::scale_color_viridis_c(name = "-log10(FDR)", option = color_palette, direction = 1) +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 15), 
                                  breaks = c(2, 5, 10, 20, 50), 
                                  labels = c("2", "5", "10", "20", "50+")) +
    ggplot2::scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50), 
                          labels = c("1", "2", "5", "10", "20", "50")) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Top ", nrow(plot_data), " significantly enriched terms (FDR < ", significance_threshold, ")"),
      x = "Fold Enrichment (log scale)",
      y = "GO Term",
      caption = paste0("Bubble size = gene count; Color = significance\n",
                      "Total terms tested: ", nrow(enrichment_results), 
                      "; Significant: ", sum(enrichment_results$p_adjusted < significance_threshold))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(size = 9, color = "gray50"),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

#' Create GO enrichment treemap
#'
#' @param enrichment_results Data frame. Output from perform_go_enrichment()
#' @param min_fold_enrichment Numeric. Minimum fold enrichment to display. Default is 2
#' @param title Character. Plot title. Default is auto-generated
#' @param color_palette Character. Color palette for significance. Default is "viridis"
#'
#' @return ggplot object
#'
#' @details
#' Creates a treemap where area represents gene count and color represents significance.
#' Useful for showing the relative importance of different GO terms.
#'
#' @examples
#' \dontrun{
#' bp_results <- perform_go_enrichment(go_data, "BP")
#' p2 <- create_go_treemap(bp_results)
#' print(p2)
#' }
#'
#' @importFrom ggplot2 ggplot aes scale_fill_viridis_c labs theme_void
#' @export
create_go_treemap <- function(enrichment_results, min_fold_enrichment = 2, title = NULL, color_palette = "viridis") {
  
  if (!requireNamespace("treemapify", quietly = TRUE)) {
    stop("treemapify package is required for this function. Install with: install.packages('treemapify')")
  }
  
  if (nrow(enrichment_results) == 0) {
    return(ggplot2::ggplot() + 
           ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No enriched terms to display", size = 5) +
           ggplot2::theme_void())
  }
  
  # Prepare data for treemap
  sig_results <- enrichment_results[enrichment_results$significance_level %in% c("significant", "highly_significant"), ]
  sig_results <- sig_results[sig_results$fold_enrichment >= min_fold_enrichment, ]
  
  treemap_data <- sig_results
  treemap_data$go_term_short <- stringr::str_trunc(treemap_data$go_term, 25)
  treemap_data$neg_log10_padj <- -log10(treemap_data$p_adjusted)
  # Add small value to avoid zero area
  treemap_data$area_value <- pmax(treemap_data$foreground_count, 1)
  
  if (nrow(treemap_data) == 0) {
    return(ggplot2::ggplot() + 
           ggplot2::annotate("text", x = 0.5, y = 0.5, 
                           label = paste0("No terms with fold enrichment >= ", min_fold_enrichment, " and FDR < 0.05"), 
                           size = 4) +
           ggplot2::theme_void())
  }
  
  # Auto-generate title if not provided
  if (is.null(title)) {
    ontology_names <- c("BP" = "Biological Process", "MF" = "Molecular Function", "CC" = "Cellular Component")
    ontology_name <- ontology_names[treemap_data$go_category[1]]
    if (is.na(ontology_name)) ontology_name <- treemap_data$go_category[1]
    title <- paste0("GO Enrichment Treemap: ", ontology_name)
  }
  
  # Create treemap
  p <- ggplot2::ggplot(treemap_data, ggplot2::aes(area = area_value, fill = neg_log10_padj, label = go_term_short)) +
    treemapify::geom_treemap(alpha = 0.8, color = "white", size = 1) +
    treemapify::geom_treemap_text(colour = "white", place = "centre", size = 10, grow = TRUE, 
                                 fontface = "bold", alpha = 0.9) +
    ggplot2::scale_fill_viridis_c(name = "-log10(FDR)", option = color_palette, direction = 1) +
    ggplot2::labs(
      title = title,
      subtitle = "Area = Gene Count; Color = Significance Level",
      caption = paste0("Showing ", nrow(treemap_data), " significantly enriched terms (FDR < 0.05, fold enrichment >= ", min_fold_enrichment, ")")
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, color = "gray50", hjust = 0.5),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Create GO enrichment summary table
#'
#' @param enrichment_results Data frame. Output from perform_go_enrichment()
#' @param max_terms Integer. Maximum number of terms to include. Default is 15
#' @param significance_filter Character vector. Significance levels to include. Default is c("significant", "highly_significant")
#'
#' @return Data frame formatted for display
#'
#' @details
#' Creates a formatted summary table of the most significant GO enrichment results,
#' suitable for inclusion in reports or manuscripts.
#'
#' @examples
#' \dontrun{
#' # Run GO enrichment workflow
#' enrich_res <- run_go_enrichment_workflow(con, candidate_vcf_file, 
#'                                          ontologies = c("BP", "MF", "CC"))
#' 
#' # Create table for specific ontology
#' bp_table <- create_go_summary_table(enrich_res$enrichment_results$BP)
#' mf_table <- create_go_summary_table(enrich_res$enrichment_results$MF)
#' 
#' # Combine all results into one table
#' all_results <- do.call(rbind, enrich_res$enrichment_results)
#' combined_table <- create_go_summary_table(all_results, max_terms = 20)
#' }
#'
#' @export
create_go_summary_table <- function(enrichment_results, max_terms = 15, 
                                   significance_filter = c("significant", "highly_significant")) {
  
  if (nrow(enrichment_results) == 0) {
    return(data.frame(Message = "No enrichment results to display"))
  }
  
  # Filter by significance
  filtered_results <- enrichment_results[enrichment_results$significance_level %in% significance_filter, ]
  
  # Arrange by p_adjusted
  filtered_results <- filtered_results[order(filtered_results$p_adjusted), ]
  
  # Take top max_terms
  if (nrow(filtered_results) > max_terms) {
    filtered_results <- filtered_results[1:max_terms, ]
  }
  
  # Create summary table
  summary_table <- data.frame(
    GO_ID = filtered_results$go_id,
    GO_Term = stringr::str_trunc(filtered_results$go_term, 60),
    Category = filtered_results$go_category,
    Genes = filtered_results$foreground_count,
    Background = filtered_results$background_count,
    Fold_Enrichment = round(filtered_results$fold_enrichment, 2),
    P_Value = format(filtered_results$p_value, scientific = TRUE, digits = 2),
    FDR = format(filtered_results$p_adjusted, scientific = TRUE, digits = 2),
    Significance = ifelse(filtered_results$significance_level == "highly_significant", "***",
                          ifelse(filtered_results$significance_level == "significant", "**",
                                 ifelse(filtered_results$significance_level == "trending", "*", ""))),
    stringsAsFactors = FALSE
  )
  
  return(summary_table)
}

#' Create multi-ontology comparison plot
#'
#' @param bp_results Data frame. Biological Process enrichment results
#' @param mf_results Data frame. Molecular Function enrichment results  
#' @param cc_results Data frame. Cellular Component enrichment results
#' @param top_n Integer. Number of top terms per category. Default is 8
#'
#' @return ggplot object showing comparison across ontologies
#'
#' @export
create_go_comparison_plot <- function(bp_results = NULL, mf_results = NULL, cc_results = NULL, top_n = 8) {
  
  # Combine results from different ontologies
  all_results <- list()
  
  if (!is.null(bp_results) && nrow(bp_results) > 0) {
    bp_filtered <- bp_results[bp_results$significance_level %in% c("significant", "highly_significant"), ]
    if (nrow(bp_filtered) > 0) {
      bp_ordered <- bp_filtered[order(bp_filtered$p_adjusted), ]
      bp_top <- if (nrow(bp_ordered) > top_n) bp_ordered[1:top_n, ] else bp_ordered
      bp_top$ontology_name <- "Biological Process"
      all_results[["BP"]] <- bp_top
    }
  }
  
  if (!is.null(mf_results) && nrow(mf_results) > 0) {
    mf_filtered <- mf_results[mf_results$significance_level %in% c("significant", "highly_significant"), ]
    if (nrow(mf_filtered) > 0) {
      mf_ordered <- mf_filtered[order(mf_filtered$p_adjusted), ]
      mf_top <- if (nrow(mf_ordered) > top_n) mf_ordered[1:top_n, ] else mf_ordered
      mf_top$ontology_name <- "Molecular Function"
      all_results[["MF"]] <- mf_top
    }
  }
  
  if (!is.null(cc_results) && nrow(cc_results) > 0) {
    cc_filtered <- cc_results[cc_results$significance_level %in% c("significant", "highly_significant"), ]
    if (nrow(cc_filtered) > 0) {
      cc_ordered <- cc_filtered[order(cc_filtered$p_adjusted), ]
      cc_top <- if (nrow(cc_ordered) > top_n) cc_ordered[1:top_n, ] else cc_ordered
      cc_top$ontology_name <- "Cellular Component"
      all_results[["CC"]] <- cc_top
    }
  }
  
  if (length(all_results) == 0) {
    return(ggplot2::ggplot() + 
           ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No significant terms in any ontology", size = 5) +
           ggplot2::theme_void())
  }
  
  # Combine all results
  combined_data <- do.call(rbind, all_results)
  combined_data$go_term_short <- stringr::str_trunc(combined_data$go_term, 45)
  combined_data$neg_log10_padj <- -log10(combined_data$p_adjusted)
  combined_data$ontology_name <- factor(combined_data$ontology_name, 
                                       levels = c("Biological Process", "Molecular Function", "Cellular Component"))
  
  # Create faceted plot
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = fold_enrichment, y = reorder(go_term_short, fold_enrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = foreground_count, color = neg_log10_padj), alpha = 0.7) +
    ggplot2::scale_color_viridis_c(name = "-log10(FDR)", option = "plasma") +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 10)) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~ontology_name, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = "GO Enrichment Across Ontologies",
      subtitle = paste0("Top ", top_n, " terms per category (FDR < 0.05)"),
      x = "Fold Enrichment (log scale)",
      y = "GO Term"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

#' Create interactive GO enrichment plot using plotly
#'
#' @param enrichment_results Data frame. Output from perform_go_enrichment()
#' @param max_terms Integer. Maximum number of terms to show. Default is 25
#' @param min_fold_enrichment Numeric. Minimum fold enrichment to display. Default is 1.5
#'
#' @return plotly object (interactive plot)
#'
#' @export
create_interactive_go_plot <- function(enrichment_results, max_terms = 25, min_fold_enrichment = 1.5) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package is required for this function. Install with: install.packages('plotly')")
  }
  
  if (nrow(enrichment_results) == 0) {
    return(NULL)
  }
  
  # Prepare data
  filtered_data <- enrichment_results[
    enrichment_results$significance_level %in% c("significant", "highly_significant") & 
    enrichment_results$fold_enrichment >= min_fold_enrichment, 
  ]
  
  if (nrow(filtered_data) == 0) {
    return(NULL)
  }
  
  ordered_data <- filtered_data[order(filtered_data$p_adjusted), ]
  plot_data <- if (nrow(ordered_data) > max_terms) ordered_data[1:max_terms, ] else ordered_data
  
  plot_data$go_term_short <- stringr::str_trunc(plot_data$go_term, 50)
  plot_data$neg_log10_padj <- -log10(plot_data$p_adjusted)
  plot_data$hover_text <- paste0(
    "<b>", plot_data$go_term, "</b><br>",
    "GO ID: ", plot_data$go_id, "<br>",
    "Genes: ", plot_data$foreground_count, "/", plot_data$total_foreground, "<br>",
    "Background: ", plot_data$background_count, "/", plot_data$total_background, "<br>",
    "Fold Enrichment: ", round(plot_data$fold_enrichment, 2), "<br>",
    "FDR: ", format(plot_data$p_adjusted, scientific = TRUE, digits = 3)
  )
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Create ggplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fold_enrichment, y = reorder(go_term_short, fold_enrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = foreground_count, color = neg_log10_padj, text = hover_text), alpha = 0.7) +
    ggplot2::scale_color_viridis_c(name = "-log10(FDR)", option = "plasma") +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 15)) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "Interactive GO Enrichment Plot",
      x = "Fold Enrichment (log scale)",
      y = "GO Term"
    ) +
    ggplot2::theme_minimal()
  
  # Convert to plotly
  interactive_plot <- plotly::ggplotly(p, tooltip = "text") %>%
    plotly::layout(
      title = list(text = "Interactive GO Enrichment Plot<br><sub>Hover for details, zoom and pan enabled</sub>"),
      showlegend = TRUE
    )
  
  return(interactive_plot)
}

#' Create KEGG pathway enrichment bubble plot
#'
#' @param enrichment_results Data frame. Output from perform_kegg_enrichment()
#' @param max_terms Integer. Maximum number of terms to show. Default is 20
#' @param min_fold_enrichment Numeric. Minimum fold enrichment to display. Default is 1.5
#' @param title Character. Plot title. Default is auto-generated
#' @param color_palette Character. Color palette for significance. Default is "plasma"
#' @param significance_threshold Numeric. FDR threshold used for significance. Default is 0.05
#'
#' @return ggplot object
#'
#' @export
create_kegg_bubble_plot <- function(enrichment_results, max_terms = 20, min_fold_enrichment = 1.5, 
                                   title = NULL, color_palette = "plasma", significance_threshold = 0.05) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("stringr package is required for this function")
  }
  
  if (nrow(enrichment_results) == 0) {
    return(ggplot2::ggplot() + 
           ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No enriched pathways to display", size = 5) +
           ggplot2::theme_void())
  }
  
  # Prepare data for plotting
  sig_results <- enrichment_results[enrichment_results$significance_level %in% c("significant", "highly_significant"), ]
  sig_results <- sig_results[sig_results$fold_enrichment >= min_fold_enrichment, ]
  sig_results <- sig_results[order(sig_results$p_adjusted), ]
  
  if (nrow(sig_results) > max_terms) {
    sig_results <- sig_results[1:max_terms, ]
  }
  
  plot_data <- sig_results
  plot_data$pathway_name_short <- stringr::str_trunc(plot_data$pathway_name, 50)
  plot_data$neg_log10_padj <- -log10(plot_data$p_adjusted)
  # Ensure minimum size for very small p-values
  plot_data$neg_log10_padj <- pmax(plot_data$neg_log10_padj, 1.3)  # -log10(0.05) = 1.3
  
  if (nrow(plot_data) == 0) {
    return(ggplot2::ggplot() + 
           ggplot2::annotate("text", x = 0.5, y = 0.5, 
                           label = paste0("No pathways with fold enrichment >= ", min_fold_enrichment, " and FDR < 0.05"), 
                           size = 4) +
           ggplot2::theme_void())
  }
  
  # Auto-generate title if not provided
  if (is.null(title)) {
    title <- "KEGG Pathway Enrichment"
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fold_enrichment, y = reorder(pathway_name_short, fold_enrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = foreground_count, color = neg_log10_padj), alpha = 0.7) +
    ggplot2::scale_color_viridis_c(name = "-log10(FDR)", option = color_palette, direction = 1) +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 15)) +
    ggplot2::scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50), 
                          labels = c("1", "2", "5", "10", "20", "50")) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Top ", nrow(plot_data), " significantly enriched pathways (FDR < ", significance_threshold, ")"),
      x = "Fold Enrichment (log scale)",
      y = "KEGG Pathway",
      caption = paste0("Bubble size = gene count; Color = significance\n",
                      "Total pathways tested: ", nrow(enrichment_results), 
                      "; Significant: ", sum(enrichment_results$p_adjusted < significance_threshold))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(size = 9, color = "gray50"),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

#' Create KEGG pathway enrichment summary table
#'
#' @param enrichment_results Data frame. Output from perform_kegg_enrichment()
#' @param max_terms Integer. Maximum number of terms to include. Default is 15
#' @param significance_filter Character vector. Significance levels to include. Default is c("significant", "highly_significant")
#'
#' @return Data frame formatted for display
#'
#' @export
create_kegg_summary_table <- function(enrichment_results, max_terms = 15, 
                                     significance_filter = c("significant", "highly_significant")) {
  
  if (nrow(enrichment_results) == 0) {
    return(data.frame(Message = "No KEGG enrichment results to display"))
  }
  
  # Filter by significance
  filtered_results <- enrichment_results[enrichment_results$significance_level %in% significance_filter, ]
  
  # Arrange by p_adjusted
  filtered_results <- filtered_results[order(filtered_results$p_adjusted), ]
  
  # Take top max_terms
  if (nrow(filtered_results) > max_terms) {
    filtered_results <- filtered_results[1:max_terms, ]
  }
  
  # Create summary table
  summary_table <- data.frame(
    Pathway_ID = filtered_results$pathway_id,
    Pathway_Name = stringr::str_trunc(filtered_results$pathway_name, 60),
    Genes = filtered_results$foreground_count,
    Background = filtered_results$background_count,
    Fold_Enrichment = round(filtered_results$fold_enrichment, 2),
    P_Value = format(filtered_results$p_value, scientific = TRUE, digits = 2),
    FDR = format(filtered_results$p_adjusted, scientific = TRUE, digits = 2),
    Significance = ifelse(filtered_results$significance_level == "highly_significant", "***",
                          ifelse(filtered_results$significance_level == "significant", "**",
                                 ifelse(filtered_results$significance_level == "trending", "*", ""))),
    stringsAsFactors = FALSE
  )
  
  return(summary_table)
}