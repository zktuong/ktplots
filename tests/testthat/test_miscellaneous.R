data(kidneyimmune)
test_that('miscellaneous works1', {
	g <- Seurat::DimPlot(kidneyimmune, color = "celltype")
	g1 <- g + small_legend() + small_guide() + small_axis() + small_grid()
	expect_true(is.ggplot(g1))
})

test_that('miscellaneous works2', {
	g <- Seurat::DimPlot(kidneyimmune, color = "celltype")
	g2 <- g + bottomleft_legend() 
	expect_true(is.ggplot(g2))
})

test_that('miscellaneous works3', {
	g <- Seurat::DimPlot(kidneyimmune, color = "celltype")
	g3 <- g + topleft_legend() 
	expect_true(is.ggplot(g3))
})

test_that('miscellaneous works4', {
	g <- Seurat::DimPlot(kidneyimmune, color = "celltype")
	g4 <- g + bottomright_legend() 
	expect_true(is.ggplot(g4))
})

test_that('miscellaneous works5', {
	g <- Seurat::DimPlot(kidneyimmune, color = "celltype")
	g5 <- g + topright_legend() 
	expect_true(is.ggplot(g5))
})
