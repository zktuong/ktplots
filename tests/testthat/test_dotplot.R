data(kidneyimmune)

test_that("geneDotPlot works", {
	p <- geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), 
		idents = "celltype", 
		split.by = 'Project', 
		standard_scale = TRUE) + theme(strip.text.x = element_text(angle=45, hjust = 0))
	expect_true(is.ggplot(p))
})