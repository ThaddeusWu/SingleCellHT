##############################################################
#  Author: [Siyuan Wu & Ulf Schmitz]                         #
#  Institution: [James Cook University]                      #
#  Date: Feb 17, 2025                                         #
#  Version: 1.2.0                                            #
##############################################################

#' Plot and analyze genes per cell distribution
#'
#' @description
#' Creates a detailed visualization of genes per cell distribution and
#' calculates suggested quality control parameters.
#'
#' @param gene_counts A matrix of gene counts (rows = genes, columns = cells)
#' @param plot_type Character, one of "hist" (histogram) or "density" (density plot)
#' @param percentile_cutoffs Numeric vector of percentiles to show as cutoff lines
#' @param return_suggestions Logical, whether to return suggested QC parameters
#'
#' @return If return_suggestions=TRUE, returns a list containing:
#'   \itemize{
#'     \item min_genes_per_cell: Suggested minimum genes per cell
#'     \item max_genes_per_cell: Suggested maximum genes per cell
#'     \item median_genes: Median number of genes per cell
#'     \item summary_stats: Additional statistics
#'   }
#'
#' @export
plot_genes_per_cell_distribution <- function(gene_counts,
                                             plot_type = "hist",
                                             percentile_cutoffs = c(0.05, 0.95),
                                             return_suggestions = TRUE) {
  # Calculate genes per cell
  n_genes <- colSums(gene_counts > 0)
  
  # Set up plotting parameters
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 4))
  
  # Create main plot
  if (plot_type == "hist") {
    hist(n_genes, 
         breaks = min(100, length(unique(n_genes))),
         main = "Distribution of Genes per Cell",
         xlab = "Number of Genes",
         ylab = "Number of Cells",
         col = "lightblue",
         border = "white")
  } else if (plot_type == "density") {
    plot(density(n_genes),
         main = "Density Distribution of Genes per Cell",
         xlab = "Number of Genes",
         ylab = "Density",
         col = "blue")
    polygon(density(n_genes), col = "lightblue", border = "blue")
  }
  
  # Add reference lines
  median_genes <- median(n_genes)
  abline(v = median_genes, col = "red", lty = 2)
  text(median_genes, par("usr")[4], 
       sprintf("Median: %d", round(median_genes)),
       pos = 3, col = "red")
  
  # Add percentile lines
  for (p in percentile_cutoffs) {
    q <- quantile(n_genes, p)
    abline(v = q, col = "blue", lty = 3)
  }
  
  # Add summary statistics
  mtext(sprintf("5th: %.1f | Mean: %.1f | Median: %.1f | SD: %.1f | CV: %.2f | 95th: %.1f",
                quantile(n_genes, 0.05), 
                mean(n_genes), 
                median(n_genes),
                sd(n_genes), 
                sd(n_genes)/mean(n_genes),
        quantile(n_genes, 0.95)), side = 3, line = 0, cex = 0.8)
  
  # Calculate suggested parameters
  if (return_suggestions) {
    suggestions <- list(
      min_genes_per_cell = as.integer(quantile(n_genes, 0.05)),
      max_genes_per_cell = as.integer(quantile(n_genes, 0.95)),
      median_genes = as.integer(median_genes),
      summary_stats = list(
        mean = mean(n_genes),
        sd = sd(n_genes),
        cv = sd(n_genes)/mean(n_genes),
        quantiles = quantile(n_genes, probs = seq(0, 1, 0.1))
      )
    )
    return(suggestions)
  }
}

#' Generate QC parameter recommendations
#'
#' @description
#' Analyzes the gene count distribution and provides detailed recommendations
#' for quality control parameter settings using different stringency levels.
#'
#' @param gene_counts A matrix of gene counts
#'
#' @return A list containing QC recommendations:
#'   \itemize{
#'     \item conservative: Conservative parameter settings
#'     \item moderate: Moderate parameter settings
#'     \item aggressive: Aggressive parameter settings
#'     \item explanation: Descriptions of each approach
#'   }
#'
#' @export
recommend_qc_parameters <- function(gene_counts) {
  # Get distribution analysis
  dist_analysis <- plot_genes_per_cell_distribution(gene_counts,
                                                    return_suggestions = TRUE)
  
  # Calculate MAD-based bounds
  median_genes <- dist_analysis$median_genes
  genes_mad <- mad(colSums(gene_counts > 0))
  conservative_min <- floor(median_genes - 3*genes_mad)
  conservative_max <- ceiling(median_genes + 3*genes_mad)
  
  # Generate recommendations
  recommendations <- list(
    conservative = list(
      min_genes_per_cell = max(100, conservative_min),  # Don't go below 100
      max_genes_per_cell = conservative_max
    ),
    moderate = list(
      min_genes_per_cell = dist_analysis$min_genes_per_cell,
      max_genes_per_cell = dist_analysis$max_genes_per_cell
    ),
    aggressive = list(
      min_genes_per_cell = as.integer(quantile(colSums(gene_counts > 0), 0.1)),
      max_genes_per_cell = as.integer(quantile(colSums(gene_counts > 0), 0.9))
    ),
    explanation = c(
      "Conservative: Uses median +/- 3 MAD, reduces risk of including poor quality cells",
      "Moderate: Uses 5th and 95th percentiles, balanced approach",
      "Aggressive: Uses 10th and 90th percentiles, retains more cells"
    )
  )
  
  return(recommendations)
}