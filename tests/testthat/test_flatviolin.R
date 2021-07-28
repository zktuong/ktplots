data(iris)
test_that('flat_violin_works', {
	p <- ggplot(iris, aes(Species, Sepal.Length)) + geom_flat_violin() + coord_flip()
	expect_true(is.ggplot(p))
})
