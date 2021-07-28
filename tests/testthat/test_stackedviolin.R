data(kidneyimmune)
features<- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

test_that('StackedVlnPlot works', {
	p <- StackedVlnPlot(kidneyimmune, features = features) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
	expect_true(is.ggplot(p))
})