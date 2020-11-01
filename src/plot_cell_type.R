#' plot_cell_categories
#' 
#' @description Produces a scatter plot of the cells in the tissue. Cells are coloured
#' categorically by phenotype. Cells not part of the phenotypes of interest will be coloured "lightgrey"
#' 
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param phenotypes_of_interest Vector of cell phenotypes to be coloured
#' @param colour_vector Vector specifying the colours of each cell phenotype
#' @import SingleCellExperiment
#' @import dplyr
#' @import ggplot2
#' @export

plot_cell_categories_CP <- function(imputed_data, phenotypes_of_interest, colour_vector) {
    
    formatted_data <- imputed_data
    
    #CHECK
    if (length(phenotypes_of_interest) != length(colour_vector)) {
        stop("The colour vector is not the same length as the phenotypes of interest")
    }
   
    
    
    formatted_data$color <- "lightgrey"
    formatted_data$Phenotype <- "OTHER"
    
    for (i in 1:length(phenotypes_of_interest)){
        
        column_name <- paste0(phenotypes_of_interest[i],"_positive")
        
        formatted_data[formatted_data[column_name]==1,]$color <- colour_vector[i]
        formatted_data[formatted_data[column_name]==1,]$Phenotype <- phenotypes_of_interest[i]
        
        # If the cells weren't found
        if (sum(formatted_data[column_name])==0) {
            print(paste(phenotype, "cells were not found"), sep="")
        }
    } 
    
    
    all_phenotypes <- c(phenotypes_of_interest, "OTHER")
    all_colours <- c(colour_vector, "lightgrey")
    
    p <- ggplot(formatted_data, aes(x = Cell.X.Position, y = Cell.Y.Position, color = Phenotype)) +
        geom_point(aes(colour = Phenotype), size = 0.5) +
        guides(alpha = F) +
        labs(colour = "Phenotypes") + 
        scale_color_manual(breaks = all_phenotypes, values=all_colours) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    
    print(p)
}