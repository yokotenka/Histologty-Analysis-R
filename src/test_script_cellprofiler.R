setwd("/Users/kentayokote/Documents/Internship/WEHI/R/Image-Analysis-Pipeline")


library(SPIAT)
source("format_csv_to_sce.R")
source("predict_positive_marker.R")
source("impute_phenotype.R")
source("pseudo_colour_environment.R")
source("tumour_deliminated_image.R")

# Path to the csv.Z
raw_data <- "Nuclei_CM_T18AH0205.csv"

##Markers used in mIHC
markers <- c("CCR7", 
             "CD103", 
             "CD3", 
             "CD4", 
             "CD45", 
             "CD45RO", 
             "CD45RA", 
             "CD49A", 
             "CD69", 
             "CD8", 
             "FOXP3", 
             "GRZB", 
             "DNA", 
             "Ki67",
             "PANCK")


# The Columns of interest in the .cpout file
##The Columns of interest in the .cpout file
intensity_columns_interest <- c("Intensity_MeanIntensity_CCR7",
                                "Intensity_MeanIntensity_CD103",
                                "Intensity_MeanIntensity_CD3",
                                "Intensity_MeanIntensity_CD4",
                                "Intensity_MeanIntensity_CD45",
                                "Intensity_MeanIntensity_CD45R0",
                                "Intensity_MeanIntensity_CD45RA",
                                "Intensity_MeanIntensity_CD49A",
                                "Intensity_MeanIntensity_CD69",
                                "Intensity_MeanIntensity_CD8",
                                "Intensity_MeanIntensity_FOXP3",
                                "Intensity_MeanIntensity_GRZB",
                                "Intensity_MeanIntensity_HEM",
                                "Intensity_MeanIntensity_KI67",
                                "Intensity_MeanIntensity_PANCK")

# Format the image to a SCE
formatted_image <- formatCsvToSce(
                                    image=raw_data,
                                    markers=markers,
                                    intensity_columns_interest=intensity_columns_interest
                                    )

# Predict whether cell is positive to marker or not
predicted <- predictPositiveMarker(
                                    formatted_image,
                                    plot_actual_cutoff = TRUE,
                                    plot_predicted_cutoff = TRUE,
                                    thresholds = NULL,
                                    tumour_marker = "PANCK",
                                    baseline_markers = c("CD45")
                                    )

# Phenotypes of interest
cellPhenotypes <- c(
    "Immune_Cell",
    "Other_Immune_Cell",
    "T_Cell",
    "Other_T_Cell",
    "CD8_T_Cell",
    "CD8_naive_T_Cell",
    "CD8_naive_T_Cell_CD69+",
    "CD8_TCM_Cell",
    "CD8_TCM_Cell_CD69+",
    "CD8_TEMRA_Cell",
    "CD8_TEMRA_Cell_CD69+",
    "CD8_Teff_Cell",
    "CD8_Teff_Cell_CD69+",
    "CD8_TEM_Cell",
    "CD8_TRM_Cell_CD69-",
    "CD8_TRM_Cell",
    "CD4_T_Cell",
    "Treg_Resting_Cell",
    "Treg_Activating_Cell",
    "Treg_Circulating_Cell",
    "Treg_TRM_Cell",
    "Treg_TRM_Cell_CD69-",
    "CD4_naive_T_Cell",
    "CD4_naive_T_Cell_CD69+",
    "CD4_TCM_Cell",
    "CD4_TCM_Cell_CD69+",
    "CD4_TEMRA_Cell",
    "CD4_TEMRA_Cell_CD69+",
    "CD4_Teff_Cell",
    "CD4_Teff_Cell_CD69+",
    "CD4_TEM_Cell",
    "CD4_TRM_Cell_CD69-",
    "CD4_TRM_Cell",
    "Epithelial_Cell",
    "Other_Cell_Type"
)

# Markers positive for the above phenotypes
positiveMarkers <- c(
    "CD45",
    "CD45",
    "CD45,CD3",
    "CD45,CD3",
    "CD45,CD3,CD8",
    "CD45,CD3,CD8,CD45RA,CCR7",
    "CD45,CD3,CD8,CD45RA,CCR7,CD69",
    "CD45,CD3,CD8,CCR7",
    "CD45,CD3,CD8,CCR7,CD69",
    "CD45,CD3,CD8,CD45RA",
    "CD45,CD3,CD8,CD45RA,CD69",
    "CD45,CD3,CD8",
    "CD45,CD3,CD8,CD69",
    "CD45,CD3,CD8,CD45RO",
    "CD45,CD3,CD8,CD45RO,CD103",
    "CD45,CD3,CD8,CD45RO,CD69",
    "CD45,CD3,CD4",
    "CD45,CD3,CD4,FOXP3,CD45RA",
    "CD45,CD3,CD4,FOXP3,CD45RO",
    "CD45,CD3,CD4,FOXP3,CD45RO",
    "CD45,CD3,CD4,FOXP3,CD45RO,CD69",
    "CD45,CD3,CD4,FOXP3,CD45RO,CD103",
    "CD45,CD3,CD4,CD45RA,CCR7",
    "CD45,CD3,CD4,CD45RA,CCR7,CD69",
    "CD45,CD3,CD4,CCR7",
    "CD45,CD3,CD4,CCR7,CD69",
    "CD45,CD3,CD4,CD45RA",
    "CD45,CD3,CD4,CD45RA,CD69",
    "CD45,CD3,CD4",
    "CD45,CD3,CD4,CD69",
    "CD45,CD3,CD4,CD45RO",
    "CD45,CD3,CD4,CD45RO,CD103",
    "CD45,CD3,CD4,CD45RO,CD69",
    "PANCK",
    "DNA"
)

# Markers negative for the above phenotype
negativeMarkers <- c(
    "PANCK",
    "PANCK,CD3",
    "PANCK",
    "PANCK,CD4,CD8",
    "PANCK,CD4,FOXP3",
    "PANCK,CD4,FOXP3,CD69,CD103",
    "PANCK,CD4,FOXP3",
    "PANCK,CD4,FOXP3,CD45RA,CD69,CD103",
    "PANCK,CD4,FOXP3,CD45RA",
    "PANCK,CD4,FOXP3,CCR7,CD69,CD103",
    "PANCK,CD4,FOXP3,CCR7",
    "PANCK,CD4,FOXP3,CD45RA,CCR7,CD45RO,CD69,CD103",
    "PANCK,CD4,FOXP3,CD45RA,CCR7,CD45RO",
    "PANCK,CD4,FOXP3,CD45RA,CCR7,CD69,CD103",
    "PANCK,CD4,FOXP3,CD45RA,CCR7,CD69",
    "PANCK,CD4,FOXP3,CD45RA,CCR7",
    "PANCK,CD8",
    "PANCK,CD8,CD45RO",
    "PANCK,CD8,CD45RA",
    "PANCK,CD8,CD45RA,CCR7,CD69,CD103",
    "PANCK,CD8,CD45RA,CCR7",
    "PANCK,CD8,CD45RA,CCR7,CD69",
    "PANCK,CD8,FOXP3,CD69,CD103",
    "PANCK,CD8,FOXP3",
    "PANCK,CD8,FOXP3,CD45RA,CD69,CD103",
    "PANCK,CD8,FOXP3,CD45RA",
    "PANCK,CD8,FOXP3,CCR7,CD69,CD103",
    "PANCK,CD8,FOXP3,CCR7",
    "PANCK,CD8,FOXP3,CCR7,CD45RA,CD45RO,CD69,CD103",
    "PANCK,CD8,FOXP3,CCR7,CD45RA,CD45RO",
    "PANCK,CD8,FOXP3,CCR7,CD45RA,CD69,CD103",
    "PANCK,CD8,FOXP3,CCR7,CD45RA,CD69",
    "PANCK,CD8,FOXP3,CCR7,CD45RA",
    "CD45,CD3,CD8,CD4,FOXP3,CD45RA,CD45RO,CCR7,CD69,CD103,CD49A",
    "CCR7,CD103,CD3,CD4,CD45,CD45RO,CD45RA,CD49A,CD69,CD8,FOXP3,GRZB,Ki67,PANCK"
)
# Impute the cell phenotypes
imputedPhenotype <- imputePhenotype(predicted$predicted_data, cellPhenotypes, positiveMarkers, negativeMarkers, markers)


phenotype_interest <- c("Immune_Cell", "Epithelial_Cell")
colour_v <- c("blue", "red")
plot_cell_categories_CP(imputedPhenotype, phenotype_interest, colour_v)



# Make the Pseudo-Colour image
# name <- "CM_T19MH0019"
# tumour_threshold <- predicted$marker_thresholds$PANCK
# 
# pseudo <- pseudo_colour_environment(
#     patient_name = name,
#     "/Volumes/T7/Processed/test2/CM_T19MH0019/V_reg_CROPPED_CM_T19MH0019_C8R2_PANCK.tif",
#     "/Volumes/T7/Processed/test2/CM_T19MH0019/NUCLEI_CROPPED_CM_T19MH0019_C8R3_HEMATOXYLIN.tif",
#     conda_env=NULL
# )
# 
# tumour_deliminated_image(pseudo_object = pseudo, 
#                          predicted$marker_thresholds$PANCK,
#                          "/Volumes/T7/Processed/test2/CM_T19MH0019/"
#                          )

