data(covid_sample_metadata)
data(covid_cpdb_meta)

file <- system.file("extdata", "covid_cpdb.tar.gz", package = "ktplots")
file.copy(file, ".")
system("tar -xzf covid_cpdb.tar.gz")

covid_sample_metadata$individual <- rep(c("A", "B", "C", "D", "E", "F"), 2)

test_that("compare_cpdb works 1", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary")
    p1 <- plot_compare_cpdb(out)
    p2 <- plot_compare_cpdb(out, cluster = TRUE, group = 'test')
    expect_true(nrow(out[[1]]) == 8339)
    expect_true(ncol(out[[1]]) == 5)
    expect_true(length(which(out[[1]]$padj < 0.05)) == 289)
    expect_true(length(which(out[[1]]$pval < 0.05)) == 941)
    expect_true(is.ggplot(p1))
    expect_true(is.ggplot(p2))
})


test_that("compare_cpdb works 2", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary", method = "t.test")
    expect_true(nrow(out[[1]]) == 8339)
    expect_true(ncol(out[[1]]) == 5)
    expect_true(length(which(out[[1]]$padj < 0.05)) == 372)
    expect_true(length(which(out[[1]]$pval < 0.05)) == 964)
})

test_that("compare_cpdb works 3", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary", formula = "~ Status_on_day_collection_summary + (1|individual)",
        method = "lmer")
    expect_true(nrow(out) == 8339)
    p <- plot_compare_cpdb(out, contrast = 'Status_on_day_collection_summarySevere', groups = 'Severe')
    expect_true(is.ggplot(p))
})


# test_that("compare_cpdb works 4", {
#     out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
#         celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
#             "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
#         groupby = "Status_on_day_collection_summary", p.adjust.mode = 'all')    
#     expect_true(length(which(out[[1]]$padj < 0.05)) == 0)
#     expect_true(length(which(out[[1]]$pval < 0.05)) == 941)
# })

# test_that("compare_cpdb works 5", {
#     covid_sample_metadata$Status_on_day_collection_summary <- c(rep('Severe', 3), rep('Healthy', 2), rep('notSevere', 2), rep('Healthy', 4), 'notSevere')
#     out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
#         celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
#             "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
#         groupby = "Status_on_day_collection_summary")
#     expect_error(plot_compare_cpdb(out))
# })
