# External
import astropy.io.fits
import os
import sys
import datetime
import glob
import time

# Internal
import image_processing
from logger import logger_obj as my_logger

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

    my_logger.info("Started main function")
    absolute_source_file_path = os.path.abspath(source_file_path)

    hdu_for_images = get_image_hdus(absolute_source_file_path, suffix_list=["b00", "d00", "f00"])

    # This is a code smell, please fix this
    bias_hdu_objects = sort_images_by_type('BIAS', hdu_for_images)
    dark_hdu_objects = sort_images_by_type('DARK', hdu_for_images)
    flat_hdu_objects = sort_images_by_type('SKYFLAT', hdu_for_images)

    bias_image_masks = image_processing.apply_bias_processing(bias_hdu_objects)

    dark_image_masks = image_processing.apply_darks_processing(dark_hdu_objects)

    flat_image_masks = image_processing.apply_flats_processing(flat_hdu_objects)

    # once you have the image masks, OR them together to get one mask
    combined_image_mask = image_processing.combine_image_masks([bias_image_masks, dark_image_masks, flat_image_masks])

    todays_date = datetime.date.today().strftime('%Y-%m-%d')
    todays_time = datetime.date.today().strftime('%H-%M-%S')

    header_dict = {
        'OBSTYPE': 'BPM',
        'DATE': todays_date,
        'TIME': todays_time
    }

    # cant directly write headers as dict object to astropy, must convert to header object first
    output_filename = os.path.join('output', "bpm-{0}-{1}.fits".format(todays_date, todays_time))
    astropy.io.fits.writeto(filename=output_filename, data=combined_image_mask, header=astropy.io.fits.Header(header_dict))

    finished_time = time.time()

    time_diff = round(finished_time - current_time, 3)
    my_logger.info("Time elapsed: {0}".format(time_diff))


def get_image_hdus(absolute_path_to_images, suffix_list, hdu_index=0,):
    """
    Takes an absolute directory that contains images and converts those images to a list of header data objects.

    :param absolute_path_to_images:
    :param hdu_to_retrieve: the index of the HDU you'd like. This will need to change to handle multi-extension fits files
    :return:
    """
    # Uses the suffix_list to construct a regex/wildcard pattern to search for files. So if your suffix list is:
    # a, b, c, this expression will be something like: '*[a|b|c].fits'
    file_wildcard_pattern = os.path.join(absolute_path_to_images, "*[{0}].fits".format('|'.join(suffix_list)))

    fits_images_list = [os.path.abspath(path) for path in glob.glob(file_wildcard_pattern)]

    hdu_list_for_images = [(astropy.io.fits.open(fits_image, mode='readonly'))[hdu_index] for fits_image in fits_images_list]

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
        if hdu.header['OBSTYPE'].strip() == image_type:
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
    my_logger.debug("Starting main script.")

    if sys.argv[1] is None:
        raise TypeError("Missing parameter -- no directory for source images was specified.")

    main(sys.argv[1])

