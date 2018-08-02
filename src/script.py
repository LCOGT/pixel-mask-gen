# External
import astropy.io.fits
import os
import sys
import datetime
import glob
import time
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


    current_time = time.time()

    logger.info("Started main function")
    absolute_source_file_path = os.path.abspath(source_file_path)

    hdu_for_images = get_image_hdus(absolute_source_file_path)

    # This is a code smell, please fix this
    bias_hdu_objects = sort_images_by_type('bias', hdu_for_images)
    dark_hdu_objects = sort_images_by_type('dark', hdu_for_images)
    flat_hdu_objects = sort_images_by_type('flat', hdu_for_images)

    bias_image_masks = image_processing.apply_bias_processing(bias_hdu_objects)

    dark_image_masks = image_processing.apply_darks_processing(dark_hdu_objects)

    flat_image_masks = image_processing.apply_flats_processing(flat_hdu_objects)

    # once  you have the image masks, OR them together to get one mask

    combined_image_mask = image_processing.combine_image_masks([bias_image_masks, dark_image_masks, flat_image_masks])

    header_dict = {
        'OBSTYPE': 'BPM',
        'DATE': datetime.date.today().strftime('%Y-%m-%d')
    }

    astropy.io.fits.writeto('final.bpm', data=combined_image_mask, header=header_dict, output_verify='exception')

    finished_time = time.time()

    time_diff = round(finished_time - current_time, 3)
    logger.info("Time elapsed: {0}".format(time_diff))




def get_image_hdus(absolute_path_to_images, hdu_to_retrieve=0):
    """
    Takes an absolute directory that contains images and converts those images to a list of header data objects.

    :param absolute_path_to_images:
    :return:
    """

    file_wildcard_pattern = os.path.join(absolute_path_to_images, "*.fits")

    fits_images_list = [os.path.abspath(path) for path in glob.glob(file_wildcard_pattern)]

    hdu_list_for_images = [(astropy.io.fits.open(fits_image, mode='readonly'))[0] for fits_image in fits_images_list]

    return hdu_list_for_images


def sort_images_by_type(image_type, hdu_list):
    """
    Takes in a list of HDU objects, and only preserves the HDU objects whose obstype match the specified image type.

    :param image_type: One of either: 'bias', 'dark', 'flat'.
    :return: A list of HDU objects
    :rtype: list
    """

    desired_hdus = []
    for hdu in hdu_list:
        if hdu.header['OBSTYPE'] == image_type:
            desired_hdus.append(hdu)

    return desired_hdus


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

