#' formatCsvToSce
#'
#' @description Formats a CellProfiler CSV into a single cell experiment. 
#' This is taken directly from format_image_to_sce from the SPIAT package and slightly modified. 
#' SPIAT: https://www.biorxiv.org/content/10.1101/2020.05.28.122614v1
#'
#' @export
#' @param image String of the path location of Cellprofiler CSV file
#' @param markers Vector containing the markers used for staining
#' @param intensity_columns_interest Vector with the names of the columns with the level of each marker.
#' Column names must match the order of the 'markers' parameter
#' @importFrom utils read.csv read.delim
#' @import SingleCellExperiment


format_csv_to_sce <- function(image, markers, intensity_columns_interest) {

    # Check extension is csv
    if (substr(image, nchar(image)-3+1, nchar(image)) != "csv"){
        stop("Please convert your file from a cpout to a csv format")
    }
    
    #read in the csv
    image <- read.csv(image)
    
    #CHECK - if image contains all the columns specified and vectors of same length
    image_colnames <- colnames(image)
    if (!all(intensity_columns_interest %in% image_colnames)) {
        stop("One or more Intensity_columns_interest not found in image")
    }
    marker_count <- length(markers)
    intensity_col_count <- length(intensity_columns_interest)
    
    if (marker_count != intensity_col_count) {
        stop("The number of columns and markers do not match")
    }
    
    # extract intensities
    intensity_of_markers <- image[,intensity_columns_interest]
    colnames(intensity_of_markers) <- markers
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    intensity_of_markers <- apply(intensity_of_markers, 2, function(x){
        as.numeric(as.character(x))
    })
    
    
    # grab relevant columns
    image <- image[,c("ObjectNumber", "Location_Center_X", "Location_Center_Y")]
    
    #rename Object.ID to Cell.ID
    colnames(image)[colnames(image) %in% c("ObjectNumber")] <- c("Cell.ID")
    
    #add "Cell_" in front of Cell.ID
    image$Cell.ID <- paste("Cell_", image$Cell.ID, sep="")
    
    colnames(image)[colnames(image) %in% c("Location_Center_X", "Location_Center_Y")] <- c("Cell.X.Position", "Cell.Y.Position")

    image <- image[,c("Cell.ID", "Cell.X.Position", "Cell.Y.Position")]


    #create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)
    
    #now create the SCE object...
    #grab the expression level, markers and cell IDs
    assay_data <- formatted_data[,markers]
    assay_rownames <- markers
    assay_colnames <- formatted_data[,"Cell.ID"]
    
    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)
    
    sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))
    
    rownames(sce) <- assay_rownames
    colnames(sce) <- assay_colnames
    
    #Assign the phenotype, X and Y positions as the colData
    coldata_Xpos <- formatted_data[,"Cell.X.Position"]
    coldata_Ypos <- formatted_data[,"Cell.Y.Position"]
    colData(sce)$Cell.X.Position <- coldata_Xpos
    colData(sce)$Cell.Y.Position <- coldata_Ypos
    
    return(sce)
}


