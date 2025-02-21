% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_scht.R
\name{.perform_cell_qc}
\alias{.perform_cell_qc}
\title{Perform cell-level quality control}
\usage{
.perform_cell_qc(gene_counts, transcript_counts, qc_params)
}
\arguments{
\item{gene_counts}{Gene-level counts matrix}

\item{transcript_counts}{Transcript-level counts matrix}

\item{qc_params}{List containing min_genes_per_cell and max_genes_per_cell thresholds}
}
\value{
A list containing:
\itemize{
\item gene_counts_filtered: Filtered gene counts matrix
\item transcript_counts_filtered: Filtered transcript counts matrix
\item n_filtered_cells: Number of cells removed
\item n_cells: Number of cells remaining
\item qc_metrics: Data frame of QC metrics per cell
}
}
\description{
Identifies and removes low-quality cells based on the number of detected
genes. Cells with too few or too many genes are considered problematic
and are filtered out.
}
\keyword{internal}
