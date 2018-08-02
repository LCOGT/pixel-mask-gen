# External
import astropy.io.fits
import numpy
import collections
import os
import yaml
import time
import datetime
import sys
import glob
import errno

import setup.LOGGER as logger

# Internal
import src.image_processing as image_processing
import src.image_object as image_object
import src.fits_utilities as fits_utilities


def main(source_file_path):
    """This function will execute when the module is run in main context.
     Loops over a range of camera IDs (0-99, must be two digits) and then tries to perform bad pixel masking on the
     corresponding camera with the specified camera prefix.

    :param filepath: the location of the config file, equivalent to sys.argv[1]: $python [arg0 = filename] \
    [arg1 = configuration file]

    :return: None.
    """

    if source_file_path == '':
        raise TypeError("Empty file path given.")


    absolute_source_file_path = os.path.abspath(source_file_path)





def get_hdulist_from_folder(absolute_path_to_images):
    """
    Takes an absolute directory that contains images and converts those images to a list of header data objects.

    :param absolute_path_to_images:
    :return:
    """

    file_wildcard_pattern = os.path.join(absolute_path_to_images, "*.fits")

    fits_images_list = [os.path.abspath(path) for path in glob.glob(file_wildcard_pattern)]

    hdu_list_for_images = [astropy.io.fits.open(fits_image, mode='readonly') for fits_image in fits_images_list]

    return hdu_list_for_images




def generate_flattened_list(list_to_flatten):
    """Takes in a doubly-nested (i.e. 2-layer deep) list and returns a flattened version where each element is a tuple

    :param list_to_flatten: A 2-layer deep list to flatten and convert to tuples
    :return: A flattened, 1D list, where each element is a tuple
    """
    flattened_list = []
    for sublist in list_to_flatten:
        for coords in sublist:
            flattened_list.append(tuple(coords))

    return flattened_list

if __name__ == '__main__':
    # Run the main script
    logger.debug("Starting main script.")

    if sys.argv[1] is None:
        raise TypeError("Missing parameter -- no directory for source images was specified.")

