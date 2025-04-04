% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_cpdb_heatmap.R
\name{plot_cpdb_heatmap}
\alias{plot_cpdb_heatmap}
\title{Plotting CellPhoneDB results as a heatmap}
\usage{
plot_cpdb_heatmap(
  pvals,
  cell_types = NULL,
  degs_analysis = FALSE,
  log1p_transform = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  border_color = "white",
  fontsize_row = 11,
  fontsize_col = 11,
  family = "Arial",
  main = "",
  treeheight_col = 0,
  treeheight_row = 0,
  low_col = "dodgerblue4",
  mid_col = "peachpuff",
  high_col = "deeppink4",
  alpha = 0.05,
  return_tables = FALSE,
  symmetrical = TRUE,
  default_sep = "\\\\|",
  ...
)
}
\arguments{
\item{pvals}{Dataframe corresponding to `pvalues.txt` or `relevant_interactions.txt` from CellPhoneDB.}

\item{cell_types}{vector of cell types to plot. If NULL, all cell types will be plotted.}

\item{degs_analysis}{Whether `CellPhoneDB` was run in `deg_analysis` mode}

\item{log1p_transform}{Whether to log1p transform the output.}

\item{show_rownames}{whether to show row names.}

\item{show_colnames}{whether to show column names.}

\item{scale}{scaling mode for pheatmap.}

\item{cluster_cols}{whether to cluster columns.}

\item{cluster_rows}{whether to cluster rows.}

\item{border_color}{border color.}

\item{fontsize_row}{row font size.}

\item{fontsize_col}{column font size.}

\item{family}{font family.}

\item{main}{plot title.}

\item{treeheight_col}{height of column dendrogram.}

\item{treeheight_row}{height of row dendrogram.}

\item{low_col}{low colour for heatmap.}

\item{mid_col}{middle colour for heatmap.}

\item{high_col}{high colour for heatmap.}

\item{alpha}{pvalue threshold to trim.}

\item{return_tables}{whether or not to return the results as a table rather than the heatmap}

\item{symmetrical}{whether or not to return as symmetrical matrix}

\item{default_sep}{the default separator used when CellPhoneDB was run.}

\item{...}{passed to pheatmap::pheatmap.}
}
\value{
pheatmap object of cellphone db output
}
\description{
Plotting CellPhoneDB results as a heatmap
}
\examples{
\donttest{
data(kidneyimmune)
data(cpdb_output2)
plot_cpdb_heatmap(pvals2)
}
}
