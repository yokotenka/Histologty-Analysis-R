#' imputePheontype

#' @description 
#' From the predicted data, this function will assign each cell with the phenotypes
#' it matches
#' 
#' @param predictedData 
#' The predicted expression of the marker from the function 
#' predict_phenotypes
#' @param phenotypeList 
#' The list of all of the phenotypes.    
#' @param positiveMarkerList 
#' The list of all markers which must be positive in the 
#' phenotypes listed in phenotypeList. Must be in the same order as phenotypeList
#' @param negativeMarkerList 
#' The list of all markers which must be negative in the 
#' phenotypes listed in phenotypeList. Must be in the same order as phenotypeList
#' @param markers 
#' The list of all markers used in the mIHC. 


impute_phenotype <- function(predictedData, phenotypeList, positiveMarkerList, negativeMarkerList, markers){
    
    # This is a dataframe of phenotypes and their regular expressions
    phenotypeRegexDf <- phenotypeRegex(phenotypeList, positiveMarkerList, negativeMarkerList, markers)
    
    # Data frame with the list of markers each cell is positive for. 
    listedMarkerData <- listPhenotype(predictedData)
    
    for (i in 1:length(phenotypeRegexDf$phenotype)){
        # Create a column with TRUE if the cell is of the particular phenotype, FALSE otherwise
        phenotypePositive <- data.frame(as.numeric(grepl(phenotypeRegexDf$regex[i], listedMarkerData$phenotype_markers)))
        colnames(phenotypePositive) <- paste0(phenotypeRegexDf$phenotype[i], "_positive")
        listedMarkerData <- cbind(listedMarkerData, phenotypePositive)
    }
    return(listedMarkerData)
} 


# Creates the regular expression for each cell phenotype given in phenotypeList
phenotypeRegex <- function(phenotypeList, positiveMarkerList, negativeMarkerList, markers){
    
    # Sorts the markers
    markers <- sort(markers)
    
    # Check if all dimensions are correct. 
    if (length(positiveMarkerList)!=length(negativeMarkerList) || 
        length(phenotypeList) !=length(negativeMarkerList) || 
        length(positiveMarkerList)!=length(phenotypeList)){
        stop("Vector sizes for the following do not agree: phenotypeList, positiveMarkerList, negativeMarkerList")
    }
    
    # Regular expression list for each phenotype
    regexList <- rep("", length(phenotypeList))
    
    for(j in 1:length(phenotypeList)){
        finalRegex <- "NA"
        for (i in 1:length(markers)){
            tempRegex <- ""
            if (markers[i] == "DNA"){
                next
            }
            
            # Positive for marker
            if (markers[i] %in% strsplit(positiveMarkerList[j],",")[[1]]){
                tempRegex <- paste(markers[i], "\\+", sep="")
            }
            # Negative for marker
            else if(markers[i] %in% strsplit(negativeMarkerList[j],",")[[1]]){
                tempRegex <- paste(markers[i], "\\-", sep="")
            }
            # Not specified
            else{
                tempRegex <- paste(markers[i], "(\\+|\\-)", sep="")
            }
            finalRegex <- paste(finalRegex, tempRegex, sep = ",")
        }
        regexList[j] <- finalRegex
    }
    
    regexList <- gsub("NA,", "", regexList)
    
    # Create a dataframe with phenotype and the associated regex. 
    phenotypeRegexPair <- data.frame(phenotype = phenotypeList, regex = regexList)
    return (phenotypeRegexPair)
}


# Lists markers from the data outputted by predict_phenotype
listPhenotype <- function(predicted_data){
    
    # Find the predicted data
    predicted_phenotype_colnames <- predicted_data[grepl("_predicted_phenotype", 
                                                         colnames(predicted_data))]
    
    # Create new column for predicted phenotype   
    predicted_data$phenotype_markers <- rep("NA", dim(predicted_data)[1])
    
    # Iterate through the predicted phenotype column and list all the phenotype markers
    # and non-phenotype markers. 
    for (predicted_marker in colnames(predicted_phenotype_colnames)){
        # The marker
        marker <- gsub("_predicted_phenotype", "", predicted_marker)
        
        # Finds all cells which are positive
        predicted_data[predicted_data[,predicted_marker] == 1,]$phenotype_markers <- 
            paste(predicted_data[predicted_data[,predicted_marker] == 1,]$phenotype_markers, 
                  paste0(marker,"+"), sep = ",") 
        
        # Finds all cells which are negative
        predicted_data[predicted_data[,predicted_marker] == 0,]$phenotype_markers <- 
            paste(predicted_data[predicted_data[,predicted_marker] == 0,]$phenotype_markers, 
                  paste0(marker,"-"), sep = ",") 
    }
    
    # Remove all the empty entries
    predicted_data$phenotype_markers <- gsub("NA,", "", predicted_data$phenotype_markers)
    predicted_data$phenotype_markers <- gsub("NA", "", predicted_data$phenotype_markers)
    return(predicted_data)
}
