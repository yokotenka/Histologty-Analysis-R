


count_cell_types <- function(patient_name, imputed_phenotype_data, phenotypes_of_interest){
    
    count <- rep(0, length(phenotypes_of_interest))
    for (i in 1:length(phenotypes_of_interest)){
        column_name <- paste0(phenotypes_of_interest[i],"_positive")
        
        count[i] <- sum(imputed_phenotype_data[column_name])
        
    }
    
    df <- data.frame(patient_name=count, row.names = phenotypes_of_interest)
}