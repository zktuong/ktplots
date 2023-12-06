#' Plotting CellPhoneDB results
#'
#' @param scdata single-cell data. can be Seurat/SingleCellExperiment object
#' @param cell_type1 Name of cell type 1. Accepts regex pattern.
#' @param cell_type2 Name of cell type 2. Accepts regex pattern.
#' @param celltype_key Column name in metadata/colData storing the celltype annotations. Values in this column should match the second column of the input `meta.txt` used for CellPhoneDB.
#' @param means Data frame corresponding to `means.txt` from CellPhoneDB.
#' @param pvals Data frame corresponding to `pvalues.txt` or `relevant_interactions.txt` from CellPhoneDB.
#' @param interaction_scores Data frame corresponding to `interaction_scores.txt` from CellPhoneDB version 5 onwards.
#' @param cellsign Data frame corresponding to `CellSign.txt` from CellPhoneDB version 5 onwards.
#' @param max_size max size of points.
#' @param keep_significant_only logical. Default is TRUE. Switch to FALSE if you want to plot all the results from cpdb.
#' @param splitby_key column name in the metadata/coldata table to split the spots by. Can only take columns with binary options. If specified, name to split by MUST be specified in the meta file provided to cpdb prior to analysis.
#' @param gene_family default = NULL. some predefined group of genes. can take one (or several) of these default options: 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche'. Also accepts name(s) of custom gene families.
#' @param custom_gene_family default = NULL. If provided, will update the gene_family function with this custom entry. Both `gene_family` (name of the custom family) and `custom_gene_family` must be specified for this to work. Provide either a data.frame with column names as name of family and genes in rows or a named likes like : list('customfamily' = c('genea', 'geneb', 'genec'))
#' @param genes default = NULL. can specify custom list of genes if gene_family is NULL
#' @param standard_scale logical. scale the expression to range from 0 to 1. NULL defaults to FALSE.
#' @param cluster_rows logical. whether or not to cluster the rows.
#' @param col_option specify plotting colours
#' @param default_style default = TRUE. Show all mean values and trace significant interactions with `higlight` colour. If FALSE, significant interactions will be presented as a white ring.
#' @param highlight_col colour for highlighting p <0.05
#' @param highlight_size stroke size for highlight if p < 0.05. if NULL, scales to -log10(pval).
#' @param max_highlight_size max size of stroke for highlight.
#' @param special_character_regex_pattern search pattern if the cell type names contains special character. NULL defaults to '/|:|\\?|\\*|\\+|[\\]|\\(|\\)'.
#' @param degs_analysis if is CellPhoneDB degs_analysis mode.
#' @param return_table whether or not to return as a table rather than to plot.
#' @param exclude_interactions if provided, the interactions will be removed from the output.
#' @param min_interaction_score Filtering the interactions shown by including only those above the given interaction score.
#' @param scale_alpha_by_interaction_scores Whether or not to filter values by the interaction score.
#' @param scale_alpha_by_cellsign Whether or not to filter the transparency of interactions by the cellsign.
#' @param filter_by_cellsign Filter out interactions with a 0 value cellsign.
#' @param ... passes arguments to grep for cell_type1 and cell_type2.
#' @return ggplot dot plot object of cellphone db output
#' @examples
#' \donttest{
#' data(kidneyimmune)
#' data(cpdb_output)
#' plot_cpdb(kidneyimmune, "B cell", "CD4T cell", "celltype", means, pvals, splitby_key = "Experiment", genes = c("CXCL13", "CD274", "CXCR5"))
#' plot_cpdb(kidneyimmune, "B cell", "CD4T cell", "celltype", means, pvals, splitby_key = "Experiment", gene_family = "chemokines")
#' }
#' @include utils.R
#' @import dplyr
#' @import viridis
#' @import ggplot2
#' @import reshape2
#' @export

plot_cpdb <- function(
    scdata, cell_type1, cell_type2, celltype_key, means, pvals,
    interaction_scores = NULL, cellsign = NULL, max_size = 8, keep_significant_only = TRUE,
    splitby_key = NULL, gene_family = NULL, custom_gene_family = NULL, genes = NULL,
    standard_scale = TRUE, cluster_rows = TRUE, col_option = viridis::viridis(50),
    default_style = TRUE, highlight_col = "red", highlight_size = NULL, max_highlight_size = 2,
    special_character_regex_pattern = NULL, degs_analysis = FALSE, return_table = FALSE,
    exclude_interactions = NULL, min_interaction_score = 0, scale_alpha_by_interaction_scores = FALSE,
    scale_alpha_by_cellsign = FALSE, filter_by_cellsign = FALSE, title = "", ...) {
  requireNamespace("SingleCellExperiment")
  requireNamespace("grDevices")
  if (is.null(special_character_regex_pattern)) {
    special_character_regex_pattern <- DEFAULT_SPEC_PAT
  }
  if (class(scdata) %in% c("SingleCellExperiment", "SummarizedExperiment")) {
    metadata <- SingleCellExperiment::colData(scdata)
  } else if (class(scdata) == "Seurat") {
    metadata <- scdata@meta.data
  }
  means_mat <- .prep_table(means)
  pvals_mat <- .prep_table(pvals)
  if (!is.null(interaction_scores)) {
    interaction_scores_mat <- .prep_table(interaction_scores)
  } else if (!is.null(cellsign)) {
    cellsign_mat <- .prep_table(cellsign)
  }
  col_start <- ifelse(colnames(pvals_mat)[DEFAULT_CLASS_COL] == "classification",
    DEFAULT_V5_COL_START, DEFAULT_COL_START
  )
  if (degs_analysis) {
    pvals_mat[, col_start:ncol(pvals_mat)] <- 1 - pvals_mat[, col_start:ncol(pvals_mat)]
  }
  # ok front load a 'dictionary' here.
  if (col_start == DEFAULT_V5_COL_START) {
    v5tmp <- reshape2::melt(means_mat, id.vars = colnames(means_mat)[1:col_start])
    row.names(v5tmp) <- paste0(v5tmp$id_cp_interaction, SPECIAL_SEP, gsub(
      "_",
      "-", v5tmp$interacting_pair
    ), SPECIAL_SEP, v5tmp$variable)
    v5tmp <- v5tmp[, c("is_integrin", "directionality", "classification")]
  }
  cell_type1 <- .sub_pattern(cell_type = cell_type1, pattern = special_character_regex_pattern)
  cell_type2 <- .sub_pattern(cell_type = cell_type2, pattern = special_character_regex_pattern)
  query_list <- .prep_query_group(
    data = means_mat, genes = genes, gene_family = gene_family,
    custom_gene_family = custom_gene_family
  )
  query <- query_list[["query"]]
  query_group <- query_list[["query_group"]]
  # prepare the cell_type query
  if (!is.null(splitby_key)) {
    labels <- paste0(metadata[[splitby_key]], "_", metadata[[celltype_key]])
    if (is.factor(metadata[[splitby_key]]) & is.factor(metadata[[celltype_key]])) {
      labels <- factor(labels, levels = paste0(
        levels(metadata[[splitby_key]]),
        "_", rep(levels(metadata[[celltype_key]]), each = length(levels(metadata[[splitby_key]])))
      ))
    } else {
      labels <- factor(labels)
    }
    labels <- levels(labels)
    groups <- factor(metadata[[splitby_key]])
    groups <- levels(groups)
    if (length(groups) > 0) {
      # the purpose for this step is to allow for special characters to be used
      # in the celltype grepping
      if (length(groups) > 1) {
        labels2 <- gsub(
          paste0(paste0(groups, "_"), collapse = "|"), "",
          labels
        )
      } else {
        labels2 <- gsub(paste0(groups, "_"), "", labels)
      }
      # this returns the indices from the labels
      ct1 <- grep(cell_type1, labels2, ...)
      ct2 <- grep(cell_type2, labels2, ...)
      c_type1 <- as.list(labels[ct1])
      c_type2 <- as.list(labels[ct2])
      c_type1 <- lapply(c_type1, .sub_pattern, pattern = special_character_regex_pattern)
      c_type2 <- lapply(c_type2, .sub_pattern, pattern = special_character_regex_pattern)
      grp <- as.list(groups)
      celltype <- list()
      for (i in 1:length(c_type1)) {
        celltype[[i]] <- .create_celltype_query(
          ctype1 = c_type1[[i]], ctype2 = c_type2,
          sep = DEFAULT_SEP
        )
        celltype[[i]] <- lapply(grp, .keep_interested_groups,
          ct = celltype[[i]],
          sep = DEFAULT_SEP
        )
      }
      for (i in 1:length(celltype)) {
        celltype[[i]] <- celltype[[i]][-which(celltype[[i]] == "")]
      }
      celltype <- lapply(celltype, unlist)
      if (any(unlist(lapply(celltype, is.null)))) {
        rm <- which(unlist(lapply(celltype, is.null)))
        celltype <- celltype[-rm]
      }
      cell_type <- do.call(paste0, list(celltype, collapse = "|"))
    } else {
      labels <- metadata[[celltype_key]]
      labels <- factor(labels)
      labels <- levels(labels)
      c_type1 <- as.list(grep(cell_type1, labels, value = TRUE, ...))
      c_type2 <- as.list(grep(cell_type2, labels, value = TRUE, ...))
      c_type1 <- lapply(c_type1, .sub_pattern, pattern = special_character_regex_pattern)
      c_type2 <- lapply(c_type2, .sub_pattern, pattern = special_character_regex_pattern)
      celltype <- list()
      for (i in 1:length(c_type1)) {
        celltype[[i]] <- .create_celltype_query(
          ctype1 = c_type1[[i]], ctype2 = c_type2,
          sep = DEFAULT_SEP
        )
      }
      cell_type <- do.call(paste0, list(celltype, collapse = "|"))
    }
  } else {
    labels <- metadata[[celltype_key]]
    labels <- factor(labels)
    labels <- levels(labels)
    c_type1 <- as.list(grep(cell_type1, labels, value = TRUE))
    c_type2 <- as.list(grep(cell_type2, labels, value = TRUE))
    c_type1 <- lapply(c_type1, .sub_pattern, pattern = special_character_regex_pattern)
    c_type2 <- lapply(c_type2, .sub_pattern, pattern = special_character_regex_pattern)
    celltype <- list()
    for (i in 1:length(c_type1)) {
      celltype[[i]] <- .create_celltype_query(
        ctype1 = c_type1[[i]], ctype2 = c_type2,
        sep = DEFAULT_SEP
      )
    }
    cell_type <- do.call(paste0, list(celltype, collapse = "|"))
  }
  if (!is.null(gene_family) & is.null(genes)) {
    if (length(gene_family) == 1) {
      means_mat <- .prep_data_querygroup_celltype1(
        .data = means_mat, .query_group = query_group,
        .gene_family = gene_family, .cell_type = cell_type, .celltype = celltype,
        ...
      )
      pvals_mat <- .prep_data_querygroup_celltype1(
        .data = pvals_mat, .query_group = query_group,
        .gene_family = gene_family, .cell_type = cell_type, .celltype = celltype,
        ...
      )
      if (!is.null(interaction_scores)) {
        interaction_scores_mat <- .prep_data_querygroup_celltype1(
          .data = interaction_scores_mat,
          .query_group = query_group, .gene_family = gene_family, .cell_type = cell_type,
          .celltype = celltype, ...
        )
      } else if (!is.null(cellsign)) {
        cellsign_mat <- .prep_data_querygroup_celltype1(
          .data = cellsign_mat,
          .query_group = query_group, .gene_family = gene_family, .cell_type = cell_type,
          .celltype = celltype, ...
        )
      }
    } else if (length(gene_family) > 1) {
      means_mat <- .prep_data_querygroup_celltype2(
        .data = means_mat, .query_group = query_group,
        .gene_family = gene_family, .cell_type = cell_type, .celltype = celltype,
        ...
      )
      pvals_mat <- .prep_data_querygroup_celltype2(
        .data = pvals_mat, .query_group = query_group,
        .gene_family = gene_family, .cell_type = cell_type, .celltype = celltype,
        ...
      )
      if (!is.null(interaction_scores)) {
        interaction_scores_mat <- .prep_data_querygroup_celltype2(
          .data = interaction_scores_mat,
          .query_group = query_group, .gene_family = gene_family, .cell_type = cell_type,
          .celltype = celltype, ...
        )
      } else if (!is.null(cellsign)) {
        cellsign_mat <- .prep_data_querygroup_celltype2(
          .data = cellsign_mat,
          .query_group = query_group, .gene_family = gene_family, .cell_type = cell_type,
          .celltype = celltype, ...
        )
      }
    }
  } else if (is.null(gene_family) & !is.null(genes) | is.null(gene_family) & is.null(genes)) {
    means_mat <- .prep_data_query_celltype(
      .data = means_mat, .query = query,
      .cell_type = cell_type, .celltype = celltype, ...
    )
    pvals_mat <- .prep_data_query_celltype(
      .data = pvals_mat, .query = query,
      .cell_type = cell_type, .celltype = celltype, ...
    )
    if (!is.null(interaction_scores)) {
      interaction_scores_mat <- .prep_data_query_celltype(
        .data = interaction_scores_mat, .query = query,
        .cell_type = cell_type, .celltype = celltype, ...
      )
    } else if (!is.null(cellsign)) {
      cellsign_mat <- cellsign_mat[, col_start:ncol(cellsign_mat)] # too difficult to code is properly?
    }
    # } else if (!is.null(cellsign)) {
    #   cellsign_mat <- .prep_data_query_celltype(
    #     .data = cellsign_mat, .query = query,
    #     .cell_type = cell_type, .celltype = celltype, ...
    #   )
    # }
  }
  if (length(means_mat) == 0) {
    stop("Please check your options for splitby_key and your celltypes.")
  } else {
    if (!all(dim(pvals_mat) == dim(means_mat))) {
      pvals_mat <- .prep_dimensions(pvals_mat, means_mat)
    }
  }
  # rearrange the columns so that it interleaves the two groups
  if (!is.null(splitby_key)) {
    if (length(groups) > 0) {
      grp <- as.list(groups)
      group_i <- lapply(grp, function(g) {
        gx <- grep(g, colnames(means_mat), ...)
        return(gx)
      })
      group_id <- do.call(c, group_i)
      means_mat <- means_mat[, as.vector(group_id), drop = FALSE]
      if (dim(pvals_mat)[2] > 0) {
        pvals_mat <- pvals_mat[, as.vector(group_id), drop = FALSE]
      } else {
        stop("No significant hits.")
      }
    }
  }
  if (keep_significant_only) {
    if (dim(pvals_mat)[2] == 0) {
      stop("No significant hits.")
    }
  }
  if (cluster_rows) {
    if (nrow(means_mat) > 2) {
      requireNamespace("stats")
      d <- stats::dist(as.data.frame(means_mat))
      h <- stats::hclust(d)
      means_mat <- means_mat[h$order, , drop = FALSE]
    }
  }
  # scaling
  if (standard_scale) {
    means_mat2 <- apply(means_mat, 1, range01)
    means_mat2 <- t(means_mat2)
  } else {
    means_mat2 <- means_mat
  }
  # remove rows that are entirely 0
  whichempty <- which(rowSums(means_mat2) == 0)
  if (length(whichempty) > 0) {
    means_mat2 <- means_mat2[whichempty, , drop = FALSE]
  }
  means_mat2 <- as.matrix(means_mat2)
  requireNamespace("reshape2")
  if (standard_scale) {
    df_means <- reshape2::melt(means_mat2, value.name = "scaled_means")
  } else {
    df_means <- reshape2::melt(means_mat2, value.name = "means")
  }
  pvals_mat2 <- as.matrix(pvals_mat)
  df_pvals <- reshape2::melt(pvals_mat2, value.name = "pvals")
  if (!is.null(interaction_scores)) {
    interaction_scores_mat2 <- as.matrix(interaction_scores_mat)
    df_interaction_scores <- reshape2::melt(interaction_scores_mat2, value.name = "interaction_scores")
  } else if (!is.null(cellsign)) {
    cellsign_mat2 <- as.matrix(cellsign_mat)
    df_cellsign <- reshape2::melt(cellsign_mat2, value.name = "cellsign")
  }
  # use dplyr left_join to combine df_means and the pvals column in df_pvals.
  # df_means and df_pvals should have the same Var1 and Var2. non-mathc should
  # fill with NA.
  df <- dplyr::left_join(df_means, df_pvals, by = c("Var1", "Var2"))
  if (!is.null(interaction_scores)) {
    df <- dplyr::left_join(df, df_interaction_scores, by = c("Var1", "Var2"))
  } else if (!is.null(cellsign)) {
    df <- dplyr::left_join(df, df_cellsign, by = c("Var1", "Var2"))
  }
  xp <- which(df$pvals == 1)
  if (length(xp) > 0) {
    df$pvals[which(df$pvals == 1)] <- NA
  }
  if (keep_significant_only) {
    # keep the entire row/ all the comparisons
    df_ <- split(df, as.character(df$Var1))
    anysig <- lapply(df_, function(x) {
      keep <- any(x$pvals < 0.05)
      return(keep)
    })
    df_ <- df_[which(unlist(anysig) == TRUE)]
    names(df_) <- NULL
    df <- do.call(rbind, df_)
  }
  df$pvals[which(df$pvals == 0)] <- 0.001
  df$pvals[which(df$pvals >= 0.05)] <- NA
  if (!is.null(splitby_key)) {
    if (length(groups) > 0) {
      grp <- as.list(groups)
      grp2 <- lapply(grp, function(i) {
        x <- paste0(i, "_")
        return(x)
      })
      searchterm <- do.call(paste, list(grp2, collapse = "|"))
      df$group <- gsub(searchterm, "", df$Var2)
    }
  } else {
    df$group <- df$Var2
  }
  if (keep_significant_only) {
    if (standard_scale) {
      if (length(df$scaled_means) == 0) {
        stop("No significant genes found and plotting will not proceed.")
      }
    } else {
      if (length(df$means) == 0) {
        stop("No significant genes found and plotting will not proceed.")
      }
    }
  }
  row.names(df) <- paste0(
    df$Var1, SPECIAL_SEP,
    df$Var2
  )
  df$Var2 <- gsub(DEFAULT_SEP, "-", df$Var2)
  df$Var1 <- gsub(paste0(".*", SPECIAL_SEP), "", df$Var1)
  final_levels <- unique(df$Var2)
  df$Var2 <- factor(df$Var2, unique(df$Var2))
  df$x_means_ <- df[, colnames(df_means)[3]]
  df$x_means_[df[, colnames(df)[4]] < 0.05] <- NA
  df$x_stroke <- df$x_means_
  df$x_stroke[!is.na(df$x_stroke)] <- 0
  df$x_stroke[is.na(df$x_stroke)] <- 2
  if (!is.null(exclude_interactions)) {
    df <- df[!df$Var1 %in% c(exclude_interactions), ]
  }
  if (!is.null(interaction_scores)) {
    df$x_means_[which(df$interaction_scores < 0)] <- NA
  } else if (!is.null(cellsign)) {
    df$cellsign[which(df$cellsign < 1)] <- 0.5
  }
  df$significant <- ifelse(df$pvals < 0.05, "yes", NA)
  if (all(is.na(df$significant))) {
    df$significant <- "no"
    highlight_col <- "#ffffff"
  }
  if (default_style) {
    df$significant[is.na(df$significant)] <- "no"
  }
  if (col_start == DEFAULT_V5_COL_START) {
    requireNamespace("tibble")
    df <- dplyr::left_join(df %>%
      tibble::rownames_to_column(), v5tmp %>%
      tibble::rownames_to_column(), by = "rowname")
    row.names(df) <- df$rowname
  }

  if (return_table) {
    return(df)
  } else {
    if (!is.null(interaction_scores)) {
      requireNamespace("dplyr")
      df <- df %>%
        dplyr::filter(interaction_scores >= min_interaction_score)
      if (scale_alpha_by_interaction_scores == TRUE) {
        if (default_style) {
          if (standard_scale) {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = scaled_means, size = scaled_means, alpha = interaction_scores
            ))
          } else {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = means, size = means, alpha = interaction_scores
            ))
          }
        } else {
          if (all(df$significant == "no")) {
            if (standard_scale) {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = scaled_means, size = scaled_means, alpha = interaction_scores
              ))
            } else {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = means, size = means, alpha = interaction_scores
              ))
            }
            default_style <- TRUE
          } else {
            highlight_col <- "#FFFFFF" # enforce this
            if (standard_scale) {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = highlight_size,
                  alpha = interaction_scores
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = x_stroke,
                  alpha = interaction_scores
                ))
              }
            } else {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = highlight_size,
                  alpha = interaction_scores
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = x_stroke, alpha = interaction_scores
                ))
              }
            }
          }
        }
      } else {
        if (default_style) {
          if (standard_scale) {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = scaled_means, size = scaled_means
            ))
          } else {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = means, size = means
            ))
          }
        } else {
          if (all(df$significant == "no")) {
            if (standard_scale) {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = scaled_means, size = scaled_means
              ))
            } else {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = means, size = means
              ))
            }
            default_style <- TRUE
          } else {
            highlight_col <- "#FFFFFF" # enforce this
            if (standard_scale) {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = highlight_size
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = x_stroke
                ))
              }
            } else {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = highlight_size
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = x_stroke
                ))
              }
            }
          }
        }
      }
    } else if (!is.null(cellsign)) {
      if (filter_by_cellsign == TRUE) {
        requireNamespace("dplyr")
        df <- df %>%
          dplyr::filter(cellsign >= 1)
      }
      if (scale_alpha_by_cellsign == TRUE) {
        if (default_style) {
          if (standard_scale) {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = scaled_means, size = scaled_means, alpha = cellsign
            ))
          } else {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = means, size = means, alpha = cellsign
            ))
          }
        } else {
          if (all(df$significant == "no")) {
            if (standard_scale) {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = scaled_means, size = scaled_means, alpha = cellsign
              ))
            } else {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = means, size = means, alpha = cellsign
              ))
            }
            default_style <- TRUE
          } else {
            highlight_col <- "#FFFFFF" # enforce this
            if (standard_scale) {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = highlight_size,
                  alpha = cellsign
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = x_stroke,
                  alpha = cellsign
                ))
              }
            } else {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = highlight_size,
                  alpha = cellsign
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = x_stroke, alpha = cellsign
                ))
              }
            }
          }
        }
      } else {
        if (default_style) {
          if (standard_scale) {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = scaled_means, size = scaled_means
            ))
          } else {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = means, size = means
            ))
          }
        } else {
          if (all(df$significant == "no")) {
            if (standard_scale) {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = scaled_means, size = scaled_means
              ))
            } else {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, color = significant,
                fill = means, size = means
              ))
            }
            default_style <- TRUE
          } else {
            highlight_col <- "#FFFFFF" # enforce this
            if (standard_scale) {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = highlight_size
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = scaled_means, size = scaled_means, stroke = x_stroke
                ))
              }
            } else {
              if (!is.null(highlight_size)) {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = highlight_size
                ))
              } else {
                g <- ggplot(df, aes(
                  x = Var2, y = Var1, fill = significant,
                  colour = means, size = means, stroke = x_stroke
                ))
              }
            }
          }
        }
      }
    } else {
      if (default_style) {
        if (standard_scale) {
          g <- ggplot(df, aes(
            x = Var2, y = Var1, color = significant, fill = scaled_means,
            size = scaled_means
          ))
        } else {
          g <- ggplot(df, aes(
            x = Var2, y = Var1, color = significant, fill = means,
            size = means
          ))
        }
      } else {
        if (all(df$significant == "no")) {
          if (standard_scale) {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = scaled_means, size = scaled_means
            ))
          } else {
            g <- ggplot(df, aes(
              x = Var2, y = Var1, color = significant,
              fill = means, size = means
            ))
          }
          default_style <- TRUE
        } else {
          highlight_col <- "#FFFFFF" # enforce this
          if (standard_scale) {
            if (!is.null(highlight_size)) {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, fill = significant,
                colour = scaled_means, size = scaled_means, stroke = highlight_size
              ))
            } else {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, fill = significant,
                colour = scaled_means, size = scaled_means, stroke = x_stroke
              ))
            }
          } else {
            if (!is.null(highlight_size)) {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, fill = significant,
                colour = means, size = means, stroke = highlight_size
              ))
            } else {
              g <- ggplot(df, aes(
                x = Var2, y = Var1, fill = significant,
                colour = means, size = means, stroke = x_stroke
              ))
            }
          }
        }
      }
    }
    g <- g + geom_point(pch = 21, na.rm = TRUE) + theme_bw() + theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 0, color = "#000000"
      ), axis.text.y = element_text(color = "#000000"),
      axis.title.x = element_blank(), axis.title.y = element_blank(), legend.direction = "vertical",
      legend.box = "horizontal"
    ) + scale_x_discrete(position = "top") + scale_radius(range = c(
      0,
      max_size
    )) + scale_linewidth(range = c(0, max_highlight_size))
    if (default_style) {
      g <- g + scale_colour_manual(
        values = c(yes = highlight_col, no = "#ffffff"),
        na.value = NA, na.translate = FALSE
      ) + guides(
        fill = guide_colourbar(
          barwidth = 4,
          label = TRUE, ticks = TRUE, draw.ulim = TRUE, draw.llim = TRUE, order = 1
        ),
        size = guide_legend(reverse = TRUE, order = 2), stroke = guide_legend(
          reverse = TRUE,
          order = 3
        )
      )
      if (length(col_option) == 1) {
        g <- g + scale_fill_gradientn(colors = (grDevices::colorRampPalette(c(
          "white",
          col_option
        )))(100), na.value = "white")
      } else {
        g <- g + scale_fill_gradientn(
          colors = c("white", (grDevices::colorRampPalette(col_option))(99)),
          na.value = "white"
        )
      }
    } else {
      g <- g + scale_fill_manual(
        values = highlight_col, na.value = "#ffffff",
        na.translate = TRUE
      ) + guides(
        colour = guide_colourbar(
          barwidth = 4,
          label = TRUE, ticks = TRUE, draw.ulim = TRUE, draw.llim = TRUE, order = 1
        ),
        size = guide_legend(reverse = TRUE, order = 2), stroke = guide_legend(
          reverse = TRUE,
          order = 3
        )
      )
      df2 <- df
      if (standard_scale) {
        df2$scaled_means[df$pvals < 0.05] <- NA
        g <- g + geom_point(aes(
          x = Var2, y = Var1, colour = scaled_means,
          size = scaled_means
        ), data = df2, inherit_aes = FALSE, na_rm = TRUE)
      } else {
        df2$means[df$pvals < 0.05] <- NA
        g <- g + geom_point(aes(x = Var2, y = Var1, colour = means, size = means),
          data = df2, inherit_aes = FALSE, na_rm = TRUE
        )
      }
      if (length(col_option) == 1) {
        g <- g + scale_colour_gradientn(colors = (grDevices::colorRampPalette(c(
          "white",
          col_option
        )))(100), na.value = "white")
      } else {
        g <- g + scale_colour_gradientn(
          colors = c("white", (grDevices::colorRampPalette(col_option))(99)),
          na.value = "white"
        )
      }
    }
    if (!is.null(interaction_scores) & (scale_alpha_by_interaction_scores ==
      TRUE)) {
      g <- g + scale_alpha_continuous(breaks = c(0, 25, 50, 75, 100))
    }
    if (!is.null(cellsign) & (scale_alpha_by_cellsign == TRUE)) {
      g <- g + scale_alpha_continuous(breaks = c(0, 1))
    }
    if (!is.null(highlight_size)) {
      g <- g + guides(stroke = "none")
    }
    if (title != "") {
      g <- g + ggtitle(title)
    } else if (!is.null(gene_family) & is.null(genes)) {
      if (length(gene_family) > 1) {
        gene_family <- paste(gene_family, collapse = ", ")
      }
      g <- g + ggtitle(gene_family)
    }
    return(g)
  }
}
