#' Plotting genes as dotplot
#' 
#' @param scdata single-cell data. can be seurat/summarizedexperiment object
#' @param idents vector holding the idents for each cell
#' @param genes genes you want to plot
#' @param split.by column name in the metadata/coldata table to split the spots by. If not provided, it will plot via idents provided.
#' @param pct.threshold float. required to keep gene expressed by minimum percentage of cells
#' @param save.plot logical. will try to save the pdf
#' @param h height of plot
#' @param w width of plot
#' @param filepath path to file, or path to folder
#' @param filename path to file
#' @param heat_cols colour gradient for the dot plot
#' @param col_limits set limits to the color gradient
#' @return ggplot dot plot object of selected genes
#' @examples
#' geneDotPlot(scdata, as.factor(scdata$split), genes = c("CD68","HLA-DRA"), save.plot = FALSE, groups = c("normal","tumor"))
#' @import dplyr
#' @import Matrix
#' @import ggplot2
#' @import reshape2
#' @export

geneDotPlot <- function(scdata, idents, genes, split.by = NULL, pct.threshold = 0.05, save.plot = FALSE, h = 5, w = 5, filepath = NULL, filename = NULL, heat_cols = rev(RColorBrewer::brewer.pal(9, "RdBu")), col_limits = NULL){
    require(ggplot2)
    require(dplyr)
    require(Matrix)
    require(reshape2)

    if (class(scdata) == "SummarizedExperiment") {
        cat("data provided is a SummarizedExperiment object", sep = "\n")
        cat("extracting expression matrix", sep = "\n")
        require(SummarizedExperiment)
        exp_mat <- SummarizedExperiment::assay(scdata)
        metadata <- SummarizedExperiment::ColData(scdata)
    } else if (class(scdata) == "Seurat") {
        cat("data provided is a Seurat object", sep = "\n")
        cat("extracting expression matrix", sep = "\n")
        exp_mat <- tryCatch(scdata@data, error = function(e) {
            tryCatch(GetAssayData(object = scdata), error = function(e) {
                stop(paste0("are you sure that your data is normalized?"))
                return(NULL)
            })
        })
        metadata <- scdata@meta.data
    }
    
    cat(paste0("attempting to subset the expression matrix to the ", length(genes), " genes provided"), sep = "\n")
    expr_mat_filtered <- exp_mat[row.names(exp_mat) %in% genes, ]

    cat(paste0("found ", dim(expr_mat_filtered)[1], " genes in the expression matrix", sep ="\n"))

    if(!is.null(split.by)){
        # cat("concatenating groups to idents", sep = "\n")
        labels = paste0(metadata[[split.by]], "_", idents)
    } else {
        cat("no groups information provided. defaulting to idents only", sep = "\n")
        labels = idents
    }

    cat("preparing the final dataframe ...", sep = "\n")
    quick_prep <- function(expr, label, groups. = NULL){
        expr <- Matrix(expr, sparse = FALSE)
        expr.df <- data.frame(label = label, t(expr), check.names = FALSE)
        
        meanExpr <- split(expr.df, expr.df$label)
        meanExpr <- lapply(meanExpr, function(x){
            x <- x[,-1]
            x <- x %>% colMeans
            return(x)
        })

        names(meanExpr) <- unique(label)
        meanExpr <- do.call(rbind, meanExpr)
        meanExpr <- scale(meanExpr)

        label.list <- as.list(unique(label))
        exp <- lapply(label.list, function(x) {
            exp_f <- expr.df %>% dplyr::filter(label == x) %>% dplyr::select(-matches("label"))
            return(exp_f)
        })

        cellNumbers <- do.call(rbind, lapply(exp, dim))[,1]
    
        pct <- list()
        pct <- lapply(exp, function(y) sapply(y, function(x) length(which(x > 0))))
        
        names(pct) <- unique(label)
        pct <- do.call(rbind, pct)
        final.pct <- pct/cellNumbers
        
        meltedMeanExpr <- melt(meanExpr)
        meltedfinal.pct <- melt(final.pct)
    
        df <- cbind(meltedMeanExpr, meltedfinal.pct$value)
        colnames(df) <- c("celltype", "gene", "scaled.mean", "pct")
    
        # add some groupings
        if(!is.null(groups.)){
            df$group <- groups.[1]
            for(i in 2:length(groups.)){
                df$group[grep(groups.[i], df$celltype)] <- groups.[i]
            }
            df$group <- factor(df$group, levels = groups.)
            remove.pattern <- paste0(groups., "_", collapse = "|")
            df$cell_type <- gsub(pattern = remove.pattern, "", df$celltype)
            df$cell_type <- as.factor(df$cell_type)
            df <- df[with(df, order(df$cell_type, df$group)), ]
        } else {
            df$cell_type <- as.factor(df$celltype)
            df$group <- as.factor(df$celltype)
            df <- df[order(df$cell_type), ]
        }
        return(df)
    }

    if(!is.null(split.by)){
        plot.df <- quick_prep(expr_mat_filtered, labels, unique(metadata[[split.by]]))
    } else {
        plot.df <- quick_prep(expr_mat_filtered, labels)
    }

    if(!is.null(pct.threshold)){
        cat(paste0("setting minimum percentage of cells expressing gene to be ", pct.threshold*100, "%"), sep ="\n")
        keep.genes <- plot.df %>% dplyr::filter(pct > pct.threshold) %>% dplyr::select(gene) %>% unique %>% unlist %>% as.character
        remove.genes <- plot.df %>% dplyr::filter(pct <= pct.threshold) %>% dplyr::select(gene) %>% unique %>% unlist %>% as.character
        cat("the following genes are removed", sep ="\n")
        print(remove.genes)
    } else {
        warning("are you sure you don't want to set a cut off?")
        keep.genes <- plot.df %>% dplyr::select(gene) %>% unique %>% unlist %>% as.character
    }

    # finally keep the genes for plotting
    plot.df.final <- plot.df[plot.df$gene %in% keep.genes, ]

    # subset the plotting objects
    doplot <- function(obj, group. = NULL, file_name = filename, file_path = filepath, dim_w, dim_h, limits. = col_limits, do.plot = save.plot){
        if(is.null(group.)){
            g <- ggplot(obj, aes(x = 0, y = gene, size = pct, colour = scaled.mean)) + 
                geom_point(pch = 16) +
                scale_y_discrete(position = "top") +
                scale_x_discrete(position = "bottom") +
                scale_colour_gradientn(colors = heat_cols, limits = limits., na.value = "grey90", oob = scales::squish) +
                scale_size_continuous(range = c(0,6), limits = c(0, 1)) +         
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),         
                    panel.border = element_blank(),
                    strip.background = element_blank()) +
                facet_grid(.~cell_type)
            } else {
            g <- ggplot(obj, aes(x = group, y = gene, size = pct, colour = scaled.mean)) + 
                geom_point(pch = 16) +
                scale_y_discrete(position = "top") +
                scale_x_discrete(position = "bottom") +
                scale_colour_gradientn(colors = heat_cols, limits = limits., na.value = "grey90", oob = scales::squish) +
                scale_size_continuous(range = c(0,6), limits = c(0, 1)) +         
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, hjust = 0),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),         
                    panel.border = element_blank(),
                    strip.background = element_blank()) +
                facet_grid(.~cell_type)
            }

        if(do.plot){
            if(is.null(file_name) && is.null(file_path)){
                out_path <- "./geneDotPlot.df"
                warning("no file name provided. saving plot to ", getwd(), "/geneDotPlot.pdf")
                ggsave("./geneDotPlot.pdf", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
            } else if(!is.null(file_name) && is.null(file_path)){
                cat(paste0("saving plot to ", file_name), sep ="\n")
                tryCatch(ggsave(file_name, plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE), error = function(e){
                    ggsave("./geneDotPlot.df", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning("file name provided is not suitable. saving as geneDotPlot.pdf")
                })
            } else if(is.null(file_name) && !is.null(file_path)){               
                cat(paste0("saving plot to ", file_path), sep ="\n")
                if(grepl('.pdf', file_path)){                   
                    ggsave(file_path, plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                } else {
                    dir.create(file_path, recursive = TRUE)
                    ggsave(paste0(file_path, "/geneDotPlot.df"), plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning(paste0("file path provided is not suitable. saving as ", file_path, "/geneDotPlot.pdf"))
                }
            } else if(!is.null(file_name) && !is.null(file_path)){
                cat(paste0("saving plot to ", paste0(file_path,"/",file_name)), sep ="\n")
                dir.create(file_path, recursive = TRUE)
                tryCatch(ggsave(paste0(file_path,"/",file_name), plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE), error = function(e){
                    ggsave("./geneDotPlot.df", plot = g, width = dim_w, height = dim_h, device = "pdf", useDingbats = FALSE)
                    warning("file path provided is not suitable. saving as geneDotPlot.pdf")
                })
            }
        }
        return(g)
    }
    gg <- doplot(plot.df.final, split.by, dim_w = w , dim_h = h)
    gg
    return(gg)
}