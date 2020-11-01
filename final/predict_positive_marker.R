#' predictPositiveMarker
#' 
#' @description  
#' This is taken directly from predict_phenotypess from the SPIAT package and slightly modified. 
#' SPIAT: https://www.biorxiv.org/content/10.1101/2020.05.28.122614v1
#' 
#' Produces a density plot showing actual and predicted cutoff of a
#' positive reading for marker expression. It also prints to the console of the
#' number of true positives (TP), true negatives (TN), false positives (FP) and
#' false negatives (FN) under the prediction. It returns a dataframe containing
#' the predicted expression status for a particular marker,
#' 
#' 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param thresholds (Optional) Vector of numbers specifying the cutoff of a positive reading.
#' The order must match the marker order, and it should be NA for DAPI.
#' @param tumour_marker String containing the tumour_marker used for the image.
#' @param baseline_markers Markers not found on tumour cells to refine the threshold
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom pracma findpeaks
#' @importFrom stats complete.cases density
#' @importFrom mmand threshold
#' @import ggplot2
#' @export

# %>% operator is in package 'magrittr' but imported by dplyr
# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

predict_positive_marker <- function(sce_object, plot_actual_cutoff = FALSE, plot_predicted_cutoff = FALSE, thresholds = NULL, tumour_marker,
                               baseline_markers) {
    
    # Markers and the associated Threshold
    marker_vector <- c()
    threshold_vector <- c()
    
    # Format from SCE to dataframe
    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column
    expression_matrix <- assay(sce_object)
    markers <- rownames(expression_matrix)
    
    # Check for tumour
    if (is.element(tumour_marker, markers) == FALSE) {
      stop("Tumour marker not found")
    }
    
    cell_ids <- colnames(expression_matrix)
    
    rownames(expression_matrix) <- NULL
    colnames(expression_matrix) <- NULL
    expression_matrix_t <- t(expression_matrix)
    expression_df <- data.frame(expression_matrix_t)
    colnames(expression_df) <- markers

    formatted_data <- cbind(formatted_data, expression_df)
    formatted_data <- formatted_data[complete.cases(formatted_data),]

    #add actual expression boolean value to formatted_data
    markers_no_tumour <- markers[markers != tumour_marker]

    for (marker in markers_no_tumour) {
        if (marker == "DNA") {
            next
        }
        #extract the marker expression column
        marker_specific_level <- formatted_data[,marker]

        #calculate the predictions
        if (!is.null(thresholds) && !is.na(thresholds[match(marker,markers)])) {
            #there is a threshold value specified for the marker, use the threshold
            marker_threshold <- thresholds[match(marker,markers)]
            print(paste("(", marker, " has threshold specified: ", as.character(marker_threshold), ")", sep=""))
            selected_valley_xcord <- NULL

            #get the threshold predictions
            predictions_by_threshold <- data.frame(threshold(marker_specific_level, level = marker_threshold))

        } else {
            #calculate the valleys
            expression_density <- density(marker_specific_level)
   
            valleys <- findpeaks(-(expression_density)$y)
            valley_ycords <- valleys[,1] * -1
            index <- match(valley_ycords, expression_density$y)
            valley_xcords <- expression_density$x[index]

            #create a df for the valley coordinates
            valley_df <- data.frame(cbind(valley_xcords, valley_ycords))

            #select the first valley that's greater than the maximum density and below 25% density
            ycord_max_density <- max(expression_density$y)
            xcord_max_density_index <- match(ycord_max_density, expression_density$y)
            xcord_max_density <- expression_density$x[xcord_max_density_index]

            density_threshold_for_valley <- 0.25 * ycord_max_density

            valley_df <- valley_df[valley_df$valley_xcords >= xcord_max_density, ]
            valley_df <- valley_df[valley_df$valley_ycords <= density_threshold_for_valley, ]

            selected_valley_xcord <- valley_df$valley_xcords[1]
            selected_valley_ycord <- valley_df$valley_ycords[1]
            
            marker_vector <- c(marker_vector, marker)
            threshold_vector<- c(threshold_vector, selected_valley_xcord)
      
            predictions_by_threshold <- data.frame(threshold(marker_specific_level, level = selected_valley_xcord))
        }
            #using the selected valley as the threshold
            
            colnames(predictions_by_threshold) <- paste(marker, "_predicted_phenotype", sep="")
            formatted_data <- cbind(formatted_data, predictions_by_threshold)
      }

    ### Prediction for tumour marker
    # Select cells that are positive for a marker that tumor cells don't have
    baseline_cells <- vector()
    
    for(marker in baseline_markers){
      temp <- formatted_data[, colnames(formatted_data) %in% c("Cell.ID", paste0(marker, "_predicted_phenotype"))]
      temp <- temp[temp[,2] == 1,]
      baseline_cells <- c(baseline_cells, temp$Cell.ID)
    }
    baseline_cells <- unique(baseline_cells)

    #Tumor marker levels in these cells

    formatted_data_baseline <- formatted_data[formatted_data$Cell.ID %in% baseline_cells,tumour_marker]
    cutoff_for_tumour <- quantile(formatted_data_baseline, 0.95)

    #extract the marker expression column
    tumour_specific_level <- formatted_data[,tumour_marker]

    #calculate the predictions
    if (!is.null(thresholds)) {
      #there is a threshold value specified for the marker, use the threshold
      marker_threshold <- thresholds[match(tumour_marker,markers)]
      print(paste("(", tumour_marker, " has threshold specified: ", as.character(marker_threshold), ")", sep=""))
      selected_valley_xcord <- NULL

      #get the threshold predictions
      predictions_by_threshold <- data.frame(threshold(tumour_specific_level, level = marker_threshold))

    } else {
      #calculate the valleys
      expression_density <- density(tumour_specific_level)
      valleys <- findpeaks(-(expression_density)$y)
      valley_ycords <- valleys[,1] * -1
      index <- match(valley_ycords, expression_density$y)
      valley_xcords <- expression_density$x[index]

      #create a df for the valley coordinates
      valley_df <- data.frame(cbind(valley_xcords, valley_ycords))
      selected_valley_xcord <- valley_df$valley_xcords[1]
      selected_valley_ycord <- valley_df$valley_ycords[1]

      #using the selected valley as the threshold if it is lower than the
      #level of expression of the tumour marker in non-tumour cells

      final_threshold <- ifelse(selected_valley_xcord < cutoff_for_tumour,
                                selected_valley_xcord, cutoff_for_tumour)
      
      marker_vector <- c(marker_vector, tumour_marker)
      threshold_vector<- c(threshold_vector, unname(final_threshold))
      
      predictions_by_threshold <- data.frame(threshold(tumour_specific_level, level = final_threshold))
      
      
    }
    colnames(predictions_by_threshold) <- paste(tumour_marker, "_predicted_phenotype", sep="")
    formatted_data <- cbind(formatted_data, predictions_by_threshold)
    predicted_columns <- formatted_data[grepl("_predicted_phenotype", colnames(formatted_data))]
    
    # Sort the columns in alphabetic order. 
    predicted_columns <- predicted_columns[, order(colnames(predicted_columns))]
    formatted_data[grepl("_predicted_phenotype", colnames(formatted_data))] <- NULL
    formatted_data <- cbind(formatted_data, predicted_columns)
    predicted_data <- formatted_data
    
    # The threshold used
    markerThresh <- data.frame(t(threshold_vector))
    colnames(markerThresh) <- marker_vector
    if (is.null(thresholds)){
      output <- list("predicted_data" = predicted_data, "marker_thresholds" = markerThresh)
      return (output)
    }
    else{
      output <- list("predicted_data" = predicted_data, "marker_thresholds" = thresholds)
      return (output)
    }
}
