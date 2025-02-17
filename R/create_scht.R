##############################################################
#  Single-Cell Hierarchical Tensor (SCHT) Creation Pipeline  #
#  Supplementary code for:                                   #
#   “A Hierarchical Representation of Isoform Expression     #
#    Reveals Transcriptomic Complexity in Single-Cells”      #
#                                                            #
#  Author: [Siyuan Wu & Ulf Schmitz]                         #
#  Institution: [James Cook University]                      #
#  Date: Feb 17, 2025                                         #
#  Version: 1.2.0                                           #
##############################################################

# This script takes raw gene- and transcript-level count matrices,
# performs a series of quality control and normalisation steps,
# and finally creates a Single-Cell Hierarchical Tensor (SCHT).

#######################
# Required Libraries  #
#######################
library(progress)


#########################
# Internal Utility Code #
#########################

#' @importFrom graphics abline hist mtext par polygon text
#' @importFrom stats density mad median na.omit quantile sd var
#' @importFrom progress progress_bar
NULL

#' Validate the SCHT input
#'
#' @description
#' Ensures input data meets required format specifications and performs
#' necessary type conversions. Validates gene counts, transcript counts,
#' and associated metadata for format consistency and completeness.
#'
#' @param gene_counts Gene-level counts as matrix or data frame
#' @param transcript_counts Transcript-level counts as matrix or data frame
#' @param transcript_info Data frame with columns: transcript_id, transcript_name, gene_id, gene_name, transcript_type gene_type
#' @param cell_info Optional data frame describing cells. Required if require_stage = TRUE
#' @param n_hvg Number of highly variable genes to select
#' @param qc_params List of quality control parameters:
#'   \itemize{
#'     \item min_genes_per_cell: Minimum number of genes per cell
#'     \item max_genes_per_cell: Maximum number of genes per cell
#'     \item min_cells_expressing: Minimum fraction of cells expressing a gene
#'     \item min_expr: Minimum expression threshold
#'   }
#' @param require_stage Logical, whether 'stage' must be present in 'cell_info'
#'
#' @return A list containing:
#'   \itemize{
#'     \item gene_counts: Validated gene counts matrix
#'     \item transcript_counts: Validated transcript counts matrix
#'   }
#'
#' @keywords internal
.validate_scht_input <- function(gene_counts,
                                 transcript_counts,
                                 transcript_info,
                                 cell_info = NULL,
                                 n_hvg,
                                 qc_params = list(
                                   # Cell QC parameters
                                   min_genes_per_cell = 200,       
                                   max_genes_per_cell = 10000,      
                                   # Gene/Transcript QC parameters
                                   min_cells_expressing = 0.02,   
                                   min_expr = 1e-6),
                                 require_stage = FALSE) {
  
  # Validate numeric parameters
  if (n_hvg < 1) stop("n_hvg must be positive")
  if (qc_params$min_genes_per_cell < 1) stop("min_genes_per_cell must be positive")
  if (qc_params$max_genes_per_cell < 1) stop("max_genes_per_cell must be positive")
  if (qc_params$min_cells_expressing <= 0 || qc_params$min_cells_expressing > 1) 
    stop("min_cells_per_gene must be between 0 and 1")
  if (qc_params$min_expr < 0) stop("min_expr must be non-negative")
  
  # Robust conversion for gene_counts
  if (is.data.frame(gene_counts)) {
    # Ensure all numeric columns
    if (!all(sapply(gene_counts, is.numeric))) {
      stop("All columns in gene_counts must be numeric")
    }
    gene_counts <- data.matrix(gene_counts)  # Using data.matrix instead of as.matrix
  }
  
  # Robust conversion for transcript_counts
  if (is.data.frame(transcript_counts)) {
    if (!all(sapply(transcript_counts, is.numeric))) {
      stop("All columns in transcript_counts must be numeric")
    }
    transcript_counts <- data.matrix(transcript_counts)
  }
  
  # Verify matrix conversion succeeded
  if (!is.matrix(gene_counts)) {
    stop("Failed to convert gene_counts to matrix format")
  }
  if (!is.matrix(transcript_counts)) {
    stop("Failed to convert transcript_counts to matrix format")
  }
  
  # Additional dimension checks
  if (length(dim(gene_counts)) != 2) {
    stop("gene_counts must be a 2-dimensional matrix")
  }
  if (length(dim(transcript_counts)) != 2) {
    stop("transcript_counts must be a 2-dimensional matrix")
  }
  
  # Check transcript_info
  if (!is.data.frame(transcript_info)) {
    stop("'transcript_info' must be a data frame.")
  }
  
  # Check cell_info if provided
  if (!is.null(cell_info)) {
    if (!is.data.frame(cell_info)) {
      stop("When provided, 'cell_info' must be a data frame.")
    }
    
    # Check for stage information
    if (require_stage && !"stage" %in% colnames(cell_info)) {
      stop("'cell_info' must contain a 'stage' column if require_stage = TRUE.")
    } else if (!require_stage && "stage" %in% colnames(cell_info)) {
      message("Found 'stage' column in cell_info but require_stage = FALSE. Stage information will not be used.")
    }
  } else {
    # cell_info is NULL
    if (require_stage) {
      stop("cell_info must be provided when require_stage = TRUE.")
    }
  }
  
  # Check transcript_info for necessary columns
  required_cols <- c("transcript_id", "transcript_name", "gene_id", "gene_name")
  if (!all(required_cols %in% colnames(transcript_info))) {
    stop("transcript_info must contain columns: ",
         paste(required_cols, collapse = ", "))
  }
  
  # Return the possibly converted matrices
  list(
    gene_counts = gene_counts,
    transcript_counts = transcript_counts
  )
  
}

#' Perform initial quality control on genes/transcripts
#'
#' @description
#' Performs initial quality control by identifying and removing low-quality
#' genes and transcripts based on expression patterns and detection rates.
#'
#' @param gene_counts Gene-level counts matrix
#' @param transcript_counts Transcript-level counts matrix
#' @param transcript_info Data frame with transcript annotations
#' @param qc_params List of quality control parameters including min_cells_expressing and min_expr
#'
#' @return A list containing:
#'   \itemize{
#'     \item gene_counts_filtered: Filtered gene counts matrix
#'     \item transcript_counts_filtered: Filtered transcript counts matrix
#'     \item transcript_info_filtered: Filtered transcript information
#'     \item n_filtered_genes: Number of genes removed
#'     \item n_filtered_transcripts: Number of transcripts removed
#'     \item qc_metrics: List of quality metrics for genes and transcripts
#'   }
#'
#' @keywords internal
.perform_initial_qc <- function(gene_counts,
                                transcript_counts,
                                transcript_info,
                                qc_params) {
  # Calculate QC metrics for genes
  n_cells_per_gene <- rowSums(gene_counts > 0)    # Number of cells expressing each gene
  mean_expr_per_gene <- rowMeans(gene_counts)     # Mean expression of each gene
  
  # Filter genes based on expression patterns
  poor_genes <- n_cells_per_gene < (qc_params$min_cells_expressing * ncol(gene_counts)) |  # Expressed in too few cells
    mean_expr_per_gene < qc_params$min_expr                                 # Too low mean expression
  
  # Filter transcripts similarly
  n_cells_per_transcript <- rowSums(transcript_counts > 0)
  mean_expr_per_transcript <- rowMeans(transcript_counts)
  
 # poor_trans <- n_cells_per_transcript < (qc_params$min_cells_expressing * ncol(transcript_counts)) |
  #  mean_expr_per_transcript < qc_params$min_expr
  
  poor_trans <- is.na(rowSums(transcript_counts > 0))
  
  # Return filtered data and QC metrics
  list(
    gene_counts_filtered = gene_counts[!poor_genes, ],
    transcript_counts_filtered = transcript_counts[!poor_trans, ],
    transcript_info_filtered = transcript_info[!poor_trans, ],
    n_filtered_genes = sum(poor_genes),
    n_filtered_transcripts = sum(poor_trans),
    qc_metrics = list(
      genes = data.frame(
        n_cells = n_cells_per_gene,
        mean_expr = mean_expr_per_gene
      ),
      transcripts = data.frame(
        n_cells = n_cells_per_transcript,
        mean_expr = mean_expr_per_transcript
      )
    )
  )
}

#' Perform cell-level quality control
#'
#' @description
#' Identifies and removes low-quality cells based on the number of detected
#' genes. Cells with too few or too many genes are considered problematic
#' and are filtered out.
#'
#' @param gene_counts Gene-level counts matrix
#' @param transcript_counts Transcript-level counts matrix
#' @param qc_params List containing min_genes_per_cell and max_genes_per_cell thresholds
#'
#' @return A list containing:
#'   \itemize{
#'     \item gene_counts_filtered: Filtered gene counts matrix
#'     \item transcript_counts_filtered: Filtered transcript counts matrix
#'     \item n_filtered_cells: Number of cells removed
#'     \item n_cells: Number of cells remaining
#'     \item qc_metrics: Data frame of QC metrics per cell
#'   }
#'
#' @keywords internal
.perform_cell_qc <- function(gene_counts,
                             transcript_counts,
                             qc_params) {
  # Calculate QC metrics
  n_genes_per_cell <- colSums(gene_counts > 0)  # Number of genes detected per cell

  # Filter cells based on multiple criteria
  poor_cells <- n_genes_per_cell < qc_params$min_genes_per_cell |  # Too few genes
    n_genes_per_cell > qc_params$max_genes_per_cell    # Too many genes (possible doublets)

  list(
    gene_counts_filtered = gene_counts[, !poor_cells, drop = FALSE],
    transcript_counts_filtered = transcript_counts[, !poor_cells, drop = FALSE],
    n_filtered_cells = sum(poor_cells),
    n_cells = ncol(gene_counts) - sum(poor_cells),
    qc_metrics = data.frame(
      n_genes = n_genes_per_cell
    )
  )
}


#' Perform CPM normalisation
#'
#' @description
#' Converts raw counts to Counts Per Million (CPM) by scaling the counts
#' in each cell to sum to one million. This normalization accounts for
#' differences in sequencing depth between cells.
#'
#' @param gene_counts_filtered Filtered gene-level counts matrix
#' @param transcript_counts_filtered Filtered transcript-level counts matrix
#'
#' @return A list containing:
#'   \itemize{
#'     \item gene_counts_norm: CPM-normalized gene counts
#'     \item transcript_counts_norm: CPM-normalized transcript counts
#'   }
#'
#' @keywords internal
#' Converts raw counts to Counts Per Million (CPM). 
.perform_normalisation <- function(gene_counts_filtered, 
                                   transcript_counts_filtered) {
  gene_counts_norm <- t(t(gene_counts_filtered) / colSums(gene_counts_filtered)) * 1e6
  transcript_counts_norm <- t(t(transcript_counts_filtered) / colSums(transcript_counts_filtered)) * 1e6
  
  list(
    gene_counts_norm = gene_counts_norm,
    transcript_counts_norm = transcript_counts_norm
  )
}


#' Perform log2 transformation
#'
#' @description
#' Applies log2(x + 1) transformation to normalized count data to reduce
#' the effect of outliers and make the data distribution more symmetric.
#'
#' @param gene_counts_norm CPM-normalized gene counts matrix
#' @param transcript_counts_norm CPM-normalized transcript counts matrix
#'
#' @return A list containing:
#'   \itemize{
#'     \item gene_mat: Log2-transformed gene expression matrix
#'     \item transcript_mat: Log2-transformed transcript expression matrix
#'   }
#'
#' @keywords internal

.perform_log_transform <- function(gene_counts_norm,
                                   transcript_counts_norm) {
  list(
    gene_mat = log2(gene_counts_norm + 1),
    transcript_mat = log2(transcript_counts_norm + 1)
  )
}


#' Select highly variable genes with progress tracking
#'
#' @description
#' Identifies highly variable genes based on their expression dispersion
#' (variance/mean ratio). Includes progress bar visualization for monitoring
#' the computation progress.
#'
#' @param gene_mat Normalized gene expression matrix
#' @param n_hvg Number of highly variable genes to select
#' @param qc_params List containing min_cells_expressing threshold
#' @param verbose Logical, whether to display progress messages
#'
#' @return A list containing:
#'   \itemize{
#'     \item var_genes: Vector of selected highly variable gene IDs
#'     \item dispersion: Named vector of dispersion values for selected genes
#'     \item mean_expr: Named vector of mean expression values for selected genes
#'   }
#'
#' @keywords internal
#' Uses variance/mean as a dispersion measure; picks top 'n_hvg' genes.
.select_highly_variable_genes <- function(gene_mat,
                                          n_hvg,
                                          qc_params,
                                          verbose = TRUE) {
  # Filter by minimum fraction of cells
  min_cells <- qc_params$min_cells_expressing * ncol(gene_mat)
  keep_genes <- rowSums(gene_mat > 0) >= min_cells
  gene_mat_filtered <- gene_mat[keep_genes, , drop = FALSE]
  
  if (verbose) {
    message("Selecting top highly variable genes (HVGs)...")
    message("Number of genes considered after min cell filter: ", nrow(gene_mat_filtered))
  }
  
  # Set up a progress bar if verbose
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "  HVG dispersion [:bar] :percent :eta",
      total = nrow(gene_mat_filtered),
      clear = FALSE, width = 60
    )
  }
  
  # Compute dispersion
  dispersion <- numeric(nrow(gene_mat_filtered))
  for (i in seq_len(nrow(gene_mat_filtered))) {
    row_vals <- gene_mat_filtered[i, ]
    dispersion[i] <- var(row_vals) / mean(row_vals)
    if (verbose) pb$tick()
  }
  
  # Rank by dispersion
  names(dispersion) <- rownames(gene_mat_filtered)
  top_genes <- names(sort(dispersion, decreasing = TRUE))[1:min(n_hvg, length(dispersion))]
  
  list(
    var_genes = top_genes,
    dispersion = dispersion[top_genes],
    mean_expr = rowMeans(gene_mat_filtered[top_genes, , drop = FALSE])
  )
}


#' Create SCHT for selected HVGs
#'
#' @description
#' Creates Single-Cell Hierarchical Tensor (SCHT) structure for the selected
#' highly variable genes by organizing their transcript-level expression data
#' into a hierarchical format.
#'
#' @param transcript_mat_final Normalized transcript expression matrix
#' @param var_genes Vector of selected highly variable gene IDs
#' @param transcript_info_filtered Filtered transcript annotation data frame
#' @param verbose Logical, whether to display progress messages
#'
#' @return A list containing:
#'   \itemize{
#'     \item SCHT: List of gene-wise isoform expression matrices
#'   }
#'
#' @keywords internals
.create_SCHT <- function(transcript_mat_final,
                                   var_genes,
                                   transcript_info_filtered,
                                   verbose = TRUE) {
  
  if (verbose) {
    message("Filtering transcripts to only those belonging to HVGs...")
    message(sprintf("Number of input var_genes: %d", length(var_genes)))
    message(sprintf("Number of unique gene names in transcript_info: %d", 
                    length(unique(transcript_info_filtered$gene_name))))
  }
  
  # Initialize empty list
  SCHT <- list()
  
  # Get unique gene IDs
  var_genes_list <- var_genes
  
  # Create progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "  Processing genes [:bar] :percent eta: :eta",
      total = length(var_genes_list),
      clear = FALSE,
      width = 60
    )
  }
  
  # Process each gene individually
  for (i in seq_along(var_genes_list)) {
    # Get transcripts for current gene
    isoform_list <- transcript_info_filtered[
      transcript_info_filtered$gene_id == var_genes_list[i],
    ]
    
    # Get gene name
    gene_name_selected <- unique(isoform_list$gene_name)
    
    # Get existing transcripts
    existing_transcripts <- isoform_list$transcript_id[
      isoform_list$transcript_id %in% rownames(transcript_mat_final)
    ]
    
    # Process if we have valid transcripts
    if (length(existing_transcripts) > 0) {
      # Extract and filter data
      result_data <- transcript_mat_final[existing_transcripts, , drop = FALSE]
      result_data <- na.omit(result_data)
      
      # Keep only transcripts with non-zero expression
      non_zero_transcripts <- result_data[rowSums(result_data) > 0, , drop = FALSE]
      
      if (nrow(non_zero_transcripts) > 0) {
        SCHT[[i]] <- non_zero_transcripts
        names(SCHT)[i] <- gene_name_selected
      }
    }
    
    # Update progress bar
    if (verbose) pb$tick()
  }
  
  # Remove any NULL elements
  SCHT <- SCHT[!sapply(SCHT, is.null)]
  
  
  # Only keep elements with more than one row (i.e., more than one isoform)
  SCHT <- SCHT[sapply(SCHT, function(x) nrow(x) > 1)]
  
  # After creating SCHT, remove columns (cells) that are all zeros
  for (gene_name in names(SCHT)) {
    non_zero_cols <- colSums(SCHT[[gene_name]]) > 0
    SCHT[[gene_name]] <- SCHT[[gene_name]][, non_zero_cols, drop = FALSE]
  }
  
  # Final status message
  if (verbose) {
    message(sprintf("Created isoform matrices for %d genes", length(SCHT)))
  }
  
  return(list(SCHT = SCHT))
}


#' Build the SCHT data structure
#'
#' @description
#' Creates the final SCHT object structure with all necessary attributes and
#' metadata. Computes and attaches basic statistics about the dataset.
#'
#' @param isoform_list List of gene-wise isoform expression matrices
#' @param transcript_info Data frame with transcript annotations
#' @param n_cells Total number of cells in the dataset
#' @param cell_info Optional data frame describing cells. Required if require_stage = TRUE
#' @param min_expr Minimum expression threshold for filtering
#' @param verbose Logical, whether to display progress messages
#'
#' @return An SCHT object with the following attributes:
#'   \itemize{
#'     \item creation_date: Timestamp of creation
#'     \item version: SCHT version number
#'     \item cell_info: Cell metadata
#'     \item stats: Basic statistics including number of genes, cells, etc.
#'   }
#'
#' @keywords internal
.build_scht_structure <- function(isoform_list,
                                  transcript_info,
                                  n_cells,
                                  cell_info = NULL,
                                  min_expr,
                                  verbose = TRUE) {
 
  scht <- isoform_list

  # Attach metadata
  attr(scht, "creation_date") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  attr(scht, "version") <- "1.0"
  attr(scht, "cell_info") <- cell_info
  
  # Compute basic statistics
  if (length(scht) > 0) {
    n_cells <- n_cells
    n_genes <- length(scht)
    total_transcripts <- sum(sapply(scht, nrow))
    mean_isoforms <- mean(sapply(scht, nrow))
    sparsity <- mean(sapply(scht, function(x) mean(x == 0)))
  } else {
    n_cells <- 0
    n_genes <- 0
    total_transcripts <- 0
    mean_isoforms <- 0
    sparsity <- NA
  }
  
  attr(scht, "stats") <- list(
    n_genes = n_genes,
    n_cells = n_cells,
    total_transcripts = total_transcripts,
    mean_isoforms = mean_isoforms,
    sparsity = sparsity
  )
  
  class(scht) <- "SCHT"
  scht
}

###################
# SCHT S3 Methods #
###################

#' Print method for SCHT objects
#'
#' @description
#' Provides a concise summary of an SCHT object, displaying key statistics
#' about genes, cells, and transcripts.
#'
#' @param x SCHT object to print
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @export
print.SCHT <- function(x, ...) {
  cat("Single-Cell Hierarchical Tensor (SCHT) object\n")
  cat("Number of genes:", attr(x, "stats")$n_genes, "\n")
  cat("Number of cells:", attr(x, "stats")$n_cells, "\n")
  cat("Total number of transcripts:", attr(x, "stats")$total_transcripts, "\n")
  cat("Mean isoforms per gene:", round(attr(x, "stats")$mean_isoforms, 2), "\n")
  cat("Overall sparsity:", round(attr(x, "stats")$sparsity * 100, 1), "%\n")
}

#' Summary method for SCHT objects
#'
#' @description
#' Generates a detailed summary of an SCHT object, including preprocessing
#' information and data characteristics.
#'
#' @param object SCHT object to summarize
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @export
summary.SCHT <- function(object, ...) {
  stats <- attr(object, "stats")
  preproc <- attr(object, "preprocessing")
  
  cat("SCHT Object Summary:\n")
  cat("--------------------\n")
  cat("  Cells:", stats$n_cells, "\n")
  cat("  Genes:", stats$n_genes, "\n")
  cat("  Total transcripts:", stats$total_transcripts, "\n")
  cat("  Mean isoforms:", round(stats$mean_isoforms, 2), "\n\n")
  
  if (!is.null(preproc)) {
    cat("Preprocessing Info:\n")
    cat("  HVGs selected:", length(preproc$hvg), "\n")
    cat("  QC-filtered genes:", preproc$qc_stats$n_filtered_genes, "\n")
    cat("  QC-filtered transcripts:", preproc$qc_stats$n_filtered_transcripts, "\n")
    cat("  QC-filtered cells:", preproc$qc_stats$n_filtered_cells, "\n\n")
  }
  
  cat("Data characteristics:\n")
  cat("  Sparsity:", round(stats$sparsity * 100, 1), "%\n")
  cat("  Created:", attr(object, "creation_date"), "\n")
}


#' Generate stage-specific SCHT structure
#'
#' @description
#' Creates stage-specific SCHT structures by partitioning the data according
#' to cell stages. Performs filtering and organization of data for each stage.
#'
#' @param scht Original SCHT object
#' @param cell_info Optional data frame describing cells. Required if require_stage = TRUE
#' @param qc_params List of quality control parameters
#'
#' @return A StageSCHT object containing:
#'   \itemize{
#'     \item stage_matrices: List of stage-specific expression matrices
#'     \item summary: Summary statistics for each stage
#'     \item data_log: Detailed log of data processing steps
#'   }
#'
#' @keywords internal
.generate_stage_scht <- function(scht, cell_info, qc_params) {
  if (!"stage" %in% colnames(cell_info)) {
    stop("cell_info must contain 'stage' column")
  }
  
  # Get unique stages
  stages <- unique(cell_info$stage)
  
  # Process each stage
  stage_scht <- list()
  data_log <- list()
  
  for (current_stage in stages) {
    # Get cells for current stage
    stage_cells <- cell_info[,1][cell_info$stage == current_stage]
    stage_matrices <- list()
    
    # Record initial stats for this stage
    initial_cells <- length(stage_cells)
    
    # Process each gene
    for (gene_name in names(scht)) {
      original_mat <- scht[[gene_name]]
      # Only select cells from current stage
      stage_mat <- original_mat[, colnames(original_mat) %in% stage_cells, drop = FALSE]
      
      # Remove cells with zero expression
      non_zero_cols <- colSums(stage_mat) > qc_params$min_expr
      if (any(non_zero_cols)) {
        stage_mat <- stage_mat[, non_zero_cols, drop = FALSE]
        stage_matrices[[gene_name]] <- stage_mat
      }
    }
    
    # Record changes
    data_log[[current_stage]] <- list(
      initial_cells = initial_cells,
      n_genes = length(stage_matrices)
    )
    
    stage_scht[[current_stage]] <- stage_matrices
  }
  
  # Create summary statistics
  summary_stats <- data.frame(
    stage = names(stage_scht),
    n_genes = sapply(stage_scht, length),
    n_cells = sapply(data_log, function(x) x$initial_cells),
    row.names = NULL
  )
  
  # Return results
  result <- list(
    stage_matrices = stage_scht,
    summary = summary_stats,
    data_log = data_log
  )
  
  class(result) <- "StageSCHT"
  return(result)
}


#' Create an integrated stage-specific SCHT structure
#'
#' @description
#' Integrates the original SCHT object with stage-specific analyses to create
#' a comprehensive data structure containing both global and stage-specific
#' information.
#'
#' @param scht_obj Original SCHT object
#' @param stage_result Results from stage-specific analysis
#'
#' @return An IntegratedSCHT object containing:
#'   \itemize{
#'     \item Original_Results: Original SCHT object
#'     \item stage_matrices: Stage-specific matrices
#'     \item Additional metadata and statistics
#'   }
#'
#' @keywords internal
.generate_integrated_result <- function(scht_obj, stage_result) {
  # Create basic structure
  final_result <- list(
    Original_Results = scht_obj,
    stage_matrices = stage_result$stage_matrices
  )
  
  # Add class and attributes
  class(final_result) <- "IntegratedSCHT"
  
  # Add metadata as attributes
  attr(final_result, "creation_date") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  attr(final_result, "stats") <- list(
    original = attr(scht_obj, "stats"),
    stage_specific = stage_result$summary
  )
  attr(final_result, "data_log") <- stage_result$data_log
  
  return(final_result)
}

##############################
# Integrated SCHT S3 Methods #
##############################

#' Print method for IntegratedSCHT objects
#'
#' @description
#' Displays a summary of an integrated SCHT object, including both original
#' and stage-specific information.
#'
#' @param x IntegratedSCHT object to print
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @export
print.IntegratedSCHT <- function(x, ...) {
  cat("Integrated SCHT Object\n")
  cat("=====================\n")
  cat("\nOriginal SCHT:\n")
  print(x$Original_Results)
  cat("\nStage-specific analysis:\n")
  cat("Number of stages:", length(x$stage_matrices), "\n")
}

#' Summary method for IntegratedSCHT objects
#'
#' @description
#' Provides a comprehensive summary of an integrated SCHT object, including
#' original SCHT summary and stage-specific analysis results.
#'
#' @param object IntegratedSCHT object to summarize
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @export
summary.IntegratedSCHT <- function(object, ...) {
  cat("Integrated SCHT Summary\n")
  cat("======================\n")
  cat("\nOriginal SCHT Summary:\n")
  cat("----------------------\n")
  summary(object$Original_Results)
  
  cat("\nStage-specific Summary:\n")
  cat("----------------------\n")
  stats <- attr(object, "stats")$stage_specific
  print(stats)
  
  cat("\nNote: The actual number of cells for each gene may vary from the shown n_cells,")
  cat("\nas cells with no expression for specific genes are removed from their respective matrices.")
  cat("\nThis cell-wise filtering is performed independently for each gene to maintain data quality")
  cat("\nand avoid spurious zero expressions in the stage-specific analyses.\n")
}



##################################
# Main SCHT Creation from Raw Data
##################################

#' Create a Single-Cell Hierarchical Tensor (SCHT) from raw count data
#'
#' @description
#' This function constructs a SCHT object from raw gene- and transcript-level
#' count matrices.
#'
#' @param gene_counts A matrix or data frame of gene-level raw counts
#' @param transcript_counts A matrix or data frame of transcript-level raw counts
#' @param transcript_info A data frame describing transcripts
#' @param cell_info Optional data frame describing cells. Required if require_stage = TRUE
#' @param n_hvg Number of highly variable genes (default: 1000)
#' @param qc_params List of quality control parameters:
#'   \itemize{
#'     \item min_genes_per_cell Minimum genes per cell
#'     \item max_genes_per_cell Maximum genes per cell
#'     \item min_cells_expressing Minimum fraction of cells expressing a gene
#'     \item min_expr Minimum expression threshold
#'   }
#' @param require_stage Logical, whether stage information is required
#' @param verbose Logical, whether to show progress messages
#'
#' @return A SCHT object
#' @export
#'
#' @examples
#' # Suppose you have raw count matrices 'gene_counts' and 'transcript_counts',
#' # along with 'transcript_info' and 'cell_info'. 
#' #
#' # scht_obj <- create_scht(
#' #   gene_counts = gene_counts,
#' #   transcript_counts = transcript_counts,
#' #   transcript_info = transcript_info,
#' #   cell_info = NULL,
#' #   n_hvg = 1000,
#' #   min_genes_per_cell = 10,
#' #   min_cells_per_gene = 0.05,
#' #   min_expr = 1e-6,
#' #   require_stage = FALSE,
#' #   verbose = TRUE
#' # )
#' #
#' # print(scht_obj)
#' # summary(scht_obj)
create_scht <- function(gene_counts,
                        transcript_counts,
                        transcript_info,
                        cell_info = NULL,
                        n_hvg = 1000,
                        qc_params = list(
                          min_genes_per_cell = 200,       
                          max_genes_per_cell = 10000,      
                          min_cells_expressing = 0.02,   
                          min_expr = 1e-6),
                        require_stage = FALSE,
                        verbose = TRUE) {
  
  # Step 0: Validate inputs
  if (verbose) message("Step 0: Validating inputs...")
  val_result <- .validate_scht_input(
    gene_counts,
    transcript_counts,
    transcript_info,
    cell_info,
    n_hvg,
    qc_params,
    require_stage
  )
  gene_counts <- val_result$gene_counts
  transcript_counts <- val_result$transcript_counts
  
  # Step 0: Input validation
  if (verbose) {
    message(sprintf("Input validation complete. 
                    Gene counts matrix: %d rows x %d cols. 
                    Transcript counts matrix: %d rows x %d cols.", 
                    nrow(val_result$gene_counts), 
                    ncol(val_result$gene_counts),
                    nrow(val_result$transcript_counts), 
                    ncol(val_result$transcript_counts)))
  }
  
  # Step 1: Initial QC
  if (verbose) message("Step 1: Performing initial quality control...")
  qc_result <- .perform_initial_qc(gene_counts, 
                                   transcript_counts, 
                                   transcript_info,
                                   qc_params)
  
  # Step 2: Cell QC
  if (verbose) message("Step 2: Filtering poor-quality cells...")
  cell_qc_result <- .perform_cell_qc(
    qc_result$gene_counts_filtered,
    qc_result$transcript_counts_filtered,
    qc_params
  )
  
  # Step 3: CPM normalisation
  if (verbose) message("Step 3: Performing CPM normalisation...")
  norm_result <- .perform_normalisation(
    cell_qc_result$gene_counts_filtered,
    cell_qc_result$transcript_counts_filtered
  )
  
  
  # Step 4: Log2 transform
  if (verbose) message("Step 4: Performing log2 transformation...")
  log_result <- .perform_log_transform(
    norm_result$gene_counts_norm,
    norm_result$transcript_counts_norm
  )
  gene_mat_final <- log_result$gene_mat
  transcript_mat_final <- log_result$transcript_mat
  
  # Step 5: Select highly variable genes
  if (verbose) message("Step 5: Selecting highly variable genes...")
  hvg_result <- .select_highly_variable_genes(
    gene_mat_final,
    n_hvg = n_hvg,
    qc_params,
    verbose = verbose
  )
  
  # Step 6: Create isoform matrix
  if (verbose) message("Step 6: Creating isoform matrix for HVGs...")
  isoform_result <- .create_SCHT(
    transcript_mat_final,
    hvg_result$var_genes,
    qc_result$transcript_info_filtered,
    verbose = verbose
  )
  
  # Step 7: Build final SCHT structure
  if (verbose) message("Step 7: Building the final SCHT structure...")
  scht_obj <- .build_scht_structure(
    isoform_result$SCHT,
    qc_result$transcript_info_filtered,
    cell_qc_result$n_cells,
    cell_info,
    qc_params,
    verbose = verbose
  )
  
  # Attach preprocessing metadata
  attr(scht_obj, "preprocessing") <- list(
    qc_stats = list(
      n_filtered_genes = qc_result$n_filtered_genes,
      n_filtered_transcripts = qc_result$n_filtered_transcripts,
      n_filtered_cells = cell_qc_result$n_filtered_cells
    ),
    hvg = hvg_result$var_genes,
    hvg_stats = list(
      n_hvg = n_hvg,
      mean_expr = hvg_result$mean_expr,
      dispersion = hvg_result$dispersion
    ),
    params = list(
      n_hvg = n_hvg,
      min_genes_per_cell = qc_params$min_genes_per_cell,
      min_cells_per_gene = qc_params$min_cells_per_gene,
      min_expr = qc_params$min_expr,
      require_stage = require_stage,
      verbose = verbose
    )
  )
  
  # Step 8: Build Stage-Specific SCHT structure (optional)
  if (require_stage && !is.null(cell_info)) {
    if (verbose) message("Step 8: Generating stage-specific SCHT structures...")
    stage_result <- .generate_stage_scht(scht_obj, 
                                         cell_info, 
                                         qc_params)
  }
  
  # Step 9: Build Integrated SCHT structure (optional)
  if (require_stage && !is.null(cell_info)) {
    if (verbose) message("Step 9: Integrating stage-specific SCHT structures...")
    scht_obj<- .generate_integrated_result(scht_obj, stage_result)
  }
  
  if (verbose) message("SCHT creation completed successfully.")
  return(scht_obj)
}

