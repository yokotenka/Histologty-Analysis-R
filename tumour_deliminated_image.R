setwd("/Users/yokote.k/Documents/Image-Analysis-Pipeline")
library(reticulate)

# In beta stage. Has many dependencies all of which can be installed by conda.
# Implementation of the class PseudoColour.

tumour_deliminated_image <- function(patient_name, tumour_path, haem_path, output_path, 
                                     tumour_threshold, conda_env=NULL, view_image=FALSE){
  
  # Source the PseudoColour.py
  source_python("PseudoColour.py")
  
  # The conda environment to be used. 
  if (conda_env){
    use_condaenv(conda_env)
  }
  
  # Instantiate a PseudoColour object
  pseudo <- PseudoColour(tumour_path, haem_path, tumour_threshold*255)
  pseudo$deliminate_tumour()
  
  # Save image to the location
  pseudo$save_all_im(patient_name, output_path)
  
  # View image 
  if (view_image){
    pseudo.plot_final_im()
  }  
}