plotCytoGenes = function (cyto_obj = cyto, numOfGenes = 10, colors = c("darkred", 
    "navyblue"), outputDir = "./") 
{
    cytoGenes <- cyto_obj$cytoGenes
    top_k <- tail(sort(cytoGenes), numOfGenes)
    bottom_k <- head(sort(cytoGenes), numOfGenes)
    final_list <- data.frame(CytoGenes = c(top_k, bottom_k))
    final_list <- cbind(final_list, Gene = rownames(final_list))
    final_list <- final_list[order(final_list$CytoGenes), ]
    p <- ggplot2::ggplot(data = final_list, ggplot2::aes(x = Gene, 
        y = CytoGenes, color = CytoGenes, fill = CytoGenes)) + 
        ggplot2::geom_bar(position = "dodge", stat = "identity") + 
        ggplot2::coord_flip() + ggplot2::scale_color_gradient(low = colors[2], 
        high = colors[1], guide = F) + ggplot2::scale_fill_gradient(low = adjustcolor(colors[2], 
        0.3), high = adjustcolor(colors[1], 0.3), guide = F) + 
        ggplot2::scale_x_discrete(limits = final_list$Gene) + 
        ggplot2::scale_y_continuous(breaks = as.numeric(formatC(round(signif(seq(round(min(cytoGenes, 
            na.rm = T) - 0.05, 1), round(max(cytoGenes, na.rm = T) + 
            0.05, 1), 0.2), 1), 1), digits = 1))) + ggplot2::ylab("Correlation with CytoTRACE") + 
        ggplot2::theme(legend.title = ggplot2::element_blank(), 
            legend.position = "none", legend.background = ggplot2::element_rect(), 
            axis.text.x = ggplot2::element_text(color = "black", 
                size = 15), axis.text.y = ggplot2::element_text(color = "black", 
                size = 15), axis.title.x = ggplot2::element_text(color = "black", 
                size = 20, margin = ggplot2::margin(t = 10, r = 0, 
                  b = 0, l = 0)), axis.title.y = ggplot2::element_text(color = "black", 
                size = 20, margin = ggplot2::margin(t = 0, r = 20, 
                  b = 0, l = 0)), axis.ticks.x = ggplot2::element_line(color = "black"), 
            axis.ticks.y = ggplot2::element_line(color = "black"), 
            axis.ticks.length = ggplot2::unit(0.2, "cm"), strip.background = ggplot2::element_blank(), 
            strip.text = ggplot2::element_text(colour = "black", 
                size = 17), axis.line = ggplot2::element_line(colour = "black"), 
            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
            panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
            plot.margin = ggplot2::margin(t = 0.5, r = 0.5, b = 0.5, 
                l = 0.5, unit = "cm"), panel.spacing.x = ggplot2::unit(1.5, 
                "lines"))
    p
return(p)
}

saveRDS(plotCytoGenes, 'Functions/plotCytoGenes.rds')
