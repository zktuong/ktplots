data(iris)
test_that('rainCloudPlot works', {
	p <- rainCloudPlot(data = iris, groupby = "Species", parameter = "Sepal.Length") + coord_flip()
	expect_true(is.ggplot(p))
})