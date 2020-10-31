
# Required packages
import cv2
import numpy as np
import scipy.ndimage 
import skimage.color
import skimage.filters
import skimage.io
import skimage.viewer
import skimage.morphology
import skimage.exposure
import skimage.filters
import skimage.measure
import skimage.util
import matplotlib.pyplot as plt



class PseudoColour:

    
    def __init__(self, tumour_path, haem_path, tumour_threshold, nuclei_path=None):
        self.tumour_path = tumour_path
        self.haem_path = haem_path
        self.nuclei_path = nuclei_path
        self.tumour_threshold = tumour_threshold


        # Nuclei Size
        self.max_nuclei_size = None

        # Original Images
        self.tumour_im = None
        self.haem_im = None

        # Mask Images
        self.tumour_mask = None
        self.stroma_mask = None
        self.background_mask = None

        # Tumour Deliminated Pseudo Coloured Image
        self.final_im = None



    ########### Utility ##########

    def adjust_colour(self, im, hue, saturation=77):
        """ Add color of the given hue to an RGB image.
        """
        hsv = cv2.cvtColor(im ,cv2.COLOR_BGR2HSV)
        hsv[:, :, 1] = saturation
        hsv[:, :, 0] = hue
        return cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)


    def apply_mask(self, im, mask):
        """ Apply mask to the RGB image
        """
        ## Need to check for rgb image
        masked = np.zeros_like(im)
        masked[:,:,0] = im[:,:,0] * mask
        masked[:,:,1] = im[:,:,1] * mask
        masked[:,:,2] = im[:,:,2] * mask
        return masked
    
    
    def enhance_contrast(self, im, lower_percentile=1, upper_percentile=99):
        # Contrast stretching
        gray_im = skimage.color.rgb2gray(im)
        p_lower = np.percentile(gray_im, lower_percentile)
        p_upper = np.percentile(gray_im, upper_percentile)
        return skimage.exposure.rescale_intensity(gray_im, in_range=(p_lower, p_upper))




    ########## Masks #############

    def create_tumour_mask(self):
        """ Create the mask for the tumour region
        """
        if (self.tumour_im is None):
            self.tumour_im = skimage.io.imread(self.tumour_path)

        tumour_bw = self.tumour_im > (self.tumour_threshold)

        if (self.max_nuclei_size is None):
            self.tumour_mask = skimage.morphology.remove_small_holes(tumour_bw, 1000)
        else:
            self.tumour_mask = skimage.morphology.remove_small_holes(tumour_bw, self.max_nuclei_size)


    def create_stroma_mask(self):
        """ Create the mask for the stroma region
        """

        if (self.background_mask is None):
            self.create_background_mask()
            
        if (self.tumour_mask is None):
            self.create_tumour_mask()

        self.stroma_mask = skimage.util.invert(np.logical_or(self.background_mask, self.tumour_mask))


    def create_background_mask(self):
        
        if (self.haem_im is None):
            self.haem_im = skimage.io.imread(self.haem_path)

        haem_contrast = self.enhance_contrast(self.haem_im)
        haem_blurred = skimage.filters.gaussian(haem_contrast, 4)

        haem_threshold = skimage.filters.threshold_otsu(haem_blurred)
        haem_bw_im = skimage.util.invert(haem_blurred > haem_threshold)

        # Find better way to filter by area
        haem_area_filtered = skimage.morphology.remove_small_objects(haem_bw_im, 1000)
        self.background_mask = skimage.util.invert(haem_area_filtered)
        


    ####### Main method
    def deliminate_tumour(self):
        """ Creates image with deliminated tumour
        """
        self.tumour_im = skimage.io.imread(self.tumour_path)
        self.haem_im = skimage.io.imread(self.haem_path)

        self.create_tumour_mask()
        self.create_background_mask()
        self.create_stroma_mask()

        # Apply the masks
        stroma_delim = self.apply_mask(self.haem_im, self.stroma_mask)
        background_delim = self.apply_mask(self.haem_im, self.background_mask)
        tumour_delim = self.apply_mask(self.haem_im, self.tumour_mask)

        
        # Change colour
        adjusted_tumour = self.adjust_colour(tumour_delim, 120)
        adjusted_stroma = self.adjust_colour(stroma_delim, 70)

        self.final_im = adjusted_tumour + adjusted_stroma + background_delim



    ########## Save image ###########
    def save_all_im(self, patient_name, output_path):
        """ Saves the tumour and stroma masks as well as the final image
        """
        # Tumour 
        tumour_fname = output_path + "/" + patient_name + "TUMOUR_REGION_MASK.tif"
        skimage.io.imsave(tumour_fname, self.tumour_mask)

        # Stroma
        stroma_fname = output_path + "/" + patient_name + "STROMA_REGION_MASK.tif"
        skimage.io.imsave(stroma_fname, self.stroma_mask)

        # Final im
        final_fname = output_path + "/" + patient_name + "DELIMINATED_TUMOUR.tif"
        skimage.io.imsave(final_fname, self.final_im)


    def plot_final_im(self):
        """ Plots the final image
        """
        skimage.io.imshow(self.final_im)
        plt.title("Deliminated Tumour Region")
        plt.show()