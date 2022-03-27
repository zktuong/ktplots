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
    expect_true(ncol(out[[1]]) == 4)
})


test_that("compare_cpdb works 2", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary", method = "wilcox")
    expect_true(nrow(out[[1]]) == 8339)
    expect_true(ncol(out[[1]]) == 4)
})

test_that("compare_cpdb works 3", {
    out <- compare_cpdb(cpdb_meta = covid_cpdb_meta, sample_metadata = covid_sample_metadata,
        celltypes = c("B_cell", "CD14", "CD16", "CD4", "CD8", "DCs", "MAIT", "NK_16hi",
            "NK_56hi", "Plasmablast", "Platelets", "Treg", "gdT", "pDC"), celltype_col = "initial_clustering",
        groupby = "Status_on_day_collection_summary", formula = "~ Status_on_day_collection_summary + (1|individual)",
        method = "lme")
    expect_true(nrow(out) == 8339)
    expect_true(ncol(out) == 7)
})