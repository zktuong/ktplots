data(kidneyimmune)
data(cpdb_output)
data(cpdb_output2)

test_that("plot_cpdb works 1",{
    p <- plot_cpdb("B cell", "CD4T cell", kidneyimmune, 'celltype', means, pvals, split.by = "Experiment", genes = c("CXCL13", "CD274", "CXCR5"))
    expect_true(is.ggplot(p))
})

test_that("plot_cpdb works 2",{
    p <- plot_cpdb("B cell", "CD4T cell", kidneyimmune, 'celltype', means, pvals, split.by = "Experiment", gene.family = 'chemokines')
    expect_true(is.ggplot(p))
})

test_that("werid characters are ok", {  
    # edit the example objects to simulate Rachel's objects
    # rename B cells to TRC+
    kidneyimmune$celltype <- gsub('B cell', 'TRC+', kidneyimmune$celltype)
    colnames(means) <- gsub('B cell', 'TRC+', colnames(means))
    colnames(pvals) <- gsub('B cell', 'TRC+', colnames(pvals))
    
    # rename NK cells to LTi-Like ILC3
    kidneyimmune$celltype <- gsub('NK cell', 'LTi-Like ILC3', kidneyimmune$celltype)
    colnames(means) <- gsub('NK cell', 'LTi-Like ILC3', colnames(means))
    colnames(pvals) <- gsub('NK cell', 'LTi-Like ILC3', colnames(pvals))
    
    # remove the original split.by tags
    colnames(means) <- gsub('Wilms2_|TxK1_|RCC1_|RCC2_|RCC3_|Wilms3_|TxK4_|VHLRCC_|Wilms1_|Teen_|Tx_|TxK3_|TxK2_|PapRCC', '', colnames(means))
    colnames(pvals) <- gsub('Wilms2_|TxK1_|RCC1_|RCC2_|RCC3_|Wilms3_|TxK4_|VHLRCC_|Wilms1_|Teen_|Tx_|TxK3_|TxK2_|PapRCC', '', colnames(pvals))
    
    # transpose and average to get rid of duplicate columns
    means_df <- as.data.frame(t(means[,c(12:ncol(means))]))
    pvals_df <- as.data.frame(t(pvals[,c(12:ncol(pvals))]))
    means_df$dup <- gsub('[.].*', '', row.names(means_df))
    pvals_df$dup <- gsub('[.].*', '', row.names(pvals_df))
    means_df <- split(means_df, means_df$dup)
    pvals_df <- split(pvals_df, pvals_df$dup)
    means_df <- t(do.call(rbind, lapply(means_df, function(x) { colMeans(x[1:ncol(x)-1])})))
    pvals_df <- t(do.call(rbind, lapply(pvals_df, function(x) { colMeans(x[1:ncol(x)-1])})))
    
    # concatenate with original first 11 columns
    newmeans <- cbind(means[,1:11], means_df)
    newpvals <- cbind(pvals[,1:11], pvals_df)
    
    # plot_cpdb
    p <- plot_cpdb(cell_type1 = 'TRC+', cell_type2 = 'LTi-Like ILC3', scdata = kidneyimmune,
        idents = "celltype", # column name where the cell ids are located in the metadata
        means = newmeans, pvals = newpvals,
        genes = c("LTB", "LTBR", "KITL", "KIT", "CCR6"))
    expect_true(is.ggplot(p))
})

test_that("plot_cpdb2 works 1",{
    sce <- Seurat::as.SingleCellExperiment(kidneyimmune)
    p <- plot_cpdb2(cell_type1 = 'B cell', cell_type2 = 'CD4T cell', 
        scdata = kidneyimmune,
        idents = 'celltype', # column name where the cell ids are located in the metadata
        means = means2, 
        pvals = pvals2,
        deconvoluted = decon2, # new options from here on specific to plot_cpdb2
        desiredInteractions = list(
            c('CD4T cell', 'B cell'), 
            c('B cell', 'CD4T cell')),
        interaction_grouping = interaction_annotation,
        edge_group_colors = c(
            "Activating" = "#e15759", 
            "Chemotaxis" = "#59a14f", 
            "Inhibitory" = "#4e79a7", 
            "Intracellular trafficking" = "#9c755f",
            "DC_development" = "#B07aa1",
            "Unknown" = NA
            ),
        node_group_colors = c(
            "CD4T cell" = "#86bc86", 
            "B cell" = "#79706e"),
        keep_significant_only = TRUE,
        standard_scale = TRUE,
        remove_self = TRUE
        )

    expect_true(is.ggplot(p))
})


test_that("plot_cpdb2 works 2",{
    sce <- Seurat::as.SingleCellExperiment(kidneyimmune)
    p <- plot_cpdb2(cell_type1 = 'B cell', cell_type2 = 'CD4T cell', 
        scdata = kidneyimmune,
        idents = 'celltype', # column name where the cell ids are located in the metadata
        split.by = 'Experiment', # column name where the grouping column is. Optional.
        means = means, 
        pvals = pvals,
        deconvoluted = decon, # new options from here on specific to plot_cpdb2
        desiredInteractions = list(
            c('CD4T cell', 'B cell'), 
            c('B cell', 'CD4T cell')),
        interaction_grouping = interaction_annotation,
        edge_group_colors = c(
            "Activating" = "#e15759", 
            "Chemotaxis" = "#59a14f", 
            "Inhibitory" = "#4e79a7", 
            "Intracellular trafficking" = "#9c755f",
            "DC_development" = "#B07aa1",
            "Unknown" = NA
            ),
        node_group_colors = c(
            "CD4T cell" = "#86bc86", 
            "B cell" = "#79706e"),
        keep_significant_only = TRUE,
        standard_scale = TRUE,
        remove_self = TRUE
        )

    expect_true(is.ggplot(p))
})