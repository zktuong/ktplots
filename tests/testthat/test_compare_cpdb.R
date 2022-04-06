data(covid_sample_metadata)
data(covid_cpdb_meta)

file <- system.file("extdata", "covid_cpdb.tar.gz", package = "ktplots")
file.copy(file, ".")
system("tar -xzf covid_cpdb.tar.gz")

covid_sample_metadata$individual <- rep(c("A", "B"), 6)

test_that("compare_cpdb works 1", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary")
    expect_true(nrow(out[[1]]) == 8339)
    expect_true(ncol(out[[1]]) == 5)
    expect_true(length(which(out[[1]]$padj < 0.05)) == 372)
    expect_true(length(which(out[[1]]$pval < 0.05)) == 964)
})


test_that("compare_cpdb works 2", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary", method = "wilcox.test")
    expect_true(nrow(out[[1]]) == 8339)
    expect_true(ncol(out[[1]]) == 5)
    expect_true(length(which(out[[1]]$padj < 0.05)) == 289)
    expect_true(length(which(out[[1]]$pval < 0.05)) == 941)
})

test_that("compare_cpdb works 3", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary", formula = "~ Status_on_day_collection_summary + (1|individual)",
        method = "lmer")
    expect_true(nrow(out) == 8339)
})