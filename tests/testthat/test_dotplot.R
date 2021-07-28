data(kidneyimmune)

test_that("geneDotPlot works", {
	p <- geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), 
		idents = "celltype", 
		split.by = 'Project', 
		standard_scale = TRUE) + theme(strip.text.x = element_text(angle=45, hjust = 0))
	expect_true(is.ggplot(p))
})

test_that("test fill works", {
	p <- geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), 
		idents = "celltype", 
		split.by = 'Project', 
		standard_scale = TRUE, 
		fill = TRUE) + theme(strip.text.x = element_text(angle=45, hjust = 0))
	expect_true(is.ggplot(p))
})

test_that("test no split", {
	p <- geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), 
		idents = "celltype", 
		standard_scale = TRUE, 
		fill = TRUE) + theme(strip.text.x = element_text(angle=45, hjust = 0))
	expect_true(is.ggplot(p))
})

test_that("test no scale1", {
	p <- geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), 
		idents = "celltype", 
		standard_scale = FALSE, 
		fill = TRUE) + theme(strip.text.x = element_text(angle=45, hjust = 0))
	expect_true(is.ggplot(p))
})

test_that("test no scale2", {
	p <- geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), 
		idents = "celltype", 
		scale = FALSE,
		standard_scale = FALSE, 
		fill = TRUE) + theme(strip.text.x = element_text(angle=45, hjust = 0))
	expect_true(is.ggplot(p))
})

test_that("test scale3", {
	p <- geneDotPlot(kidneyimmune, genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), 
		idents = "celltype") + theme(strip.text.x = element_text(angle=45, hjust = 0))
	expect_true(is.ggplot(p))
})