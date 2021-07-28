data(kidneyimmune)
data(cpdb_output)
test_that("geneDotPlot works 1",{
	p <- plot_cpdb("B cell", "CD4T cell", kidneyimmune, 'celltype', means, pvals, split.by = "Experiment", genes = c("CXCL13", "CD274", "CXCR5"))
	expect_true(is.ggplot(p))
})

test_that("geneDotPlot works 2",{
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
