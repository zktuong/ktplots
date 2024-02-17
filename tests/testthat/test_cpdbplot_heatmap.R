data(kidneyimmune)
data(cpdb_output)
data(cpdb_output2)
data(cpdb_output_v5)


test_that("plot_cpdb_heatmap works", {
    p <- plot_cpdb_heatmap(pvals2)
    expect_that(class(p), equals("pheatmap"))
})

test_that("plot_cpdb_heatmap return table", {
    p <- plot_cpdb_heatmap(pvals2, return_tables = TRUE)
    expect_true(is.list(p))
    expect_true(is.matrix(p[[1]]))
    expect_true(is.data.frame(p[[2]]))
})

test_that("plot_cpdb_heatmap works for log", {
    p <- plot_cpdb_heatmap(pvals2, log1p_transform = TRUE)
    expect_that(class(p), equals("pheatmap"))
})

test_that("plot_cpdb_heatmap return table", {
    p <- plot_cpdb_heatmap(pvals2, return_tables = TRUE, log1p_transform = TRUE)
    expect_true(is.list(p))
    expect_true(is.matrix(p[[1]]))
    expect_true(is.data.frame(p[[2]]))
})

test_that("plot_cpdb_heatmap works2", {
    p <- plot_cpdb_heatmap(relevant_interactions_v5, cell_types = c(
        "iEVT", "PV MYH11",
        "PV STEAP4", "EVT_1"
    ), degs_analysis = TRUE)
    expect_that(class(p), equals("pheatmap"))
})

test_that("plot_cpdb_heatmap works2", {
    p <- plot_cpdb_heatmap(relevant_interactions_v5, cell_types = c(
        "iEVT", "PV MYH11",
        "PV STEAP4", "EVT_1"
    ), degs_analysis = FALSE)
    expect_that(class(p), equals("pheatmap"))
})


test_that("plot_cpdb_heatmap return table not symmetrical", {
    p <- plot_cpdb_heatmap(relevant_interactions_v5, cell_types = c(
        "iEVT", "PV MYH11",
        "PV STEAP4", "EVT_1"
    ), degs_analysis = TRUE, return_tables = TRUE, symmetrical = FALSE)
    expect_true(is.list(p))
    expect_true(is.matrix(p[[1]]))
    expect_true(is.data.frame(p[[2]]))
})
