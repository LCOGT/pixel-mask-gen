# External
import astropy.io.fits
import os
import sys
import datetime
import glob
import time
import logging

# Internal
import src.image_processing as image_processing
import src.logger as logger

def main(source_file_path):
    r"""This function will only execute when the module is run in main context.

    This function does the following:

    Read in the images in the directory given via the command line parameter.

    Depending on the image type, apply the correct type of processing to it.

    Combine the masks for each type into one mask.

    Write that to a FITS file.

    :param source_file_path: The 2nd command-line parameter, which is an absolute path to the root/top-level directory\
    of where the images will be stored.

    :rtype: None

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

    today = datetime.datetime.now()
    format_strings = ["%Y-%m-%d", "-%H%M%S"]
    today_date, today_time = tuple([today.strftime(str_format) for str_format in format_strings])

    header_dict = {
        'OBSTYPE': 'BPM',
        'DATE': today.strftime("%Y-%m-%d"),
        'TIME': today.strftime("-%H%M%S")
    }

    # cant directly write headers as dict object to astropy, must convert to header object first
    output_filename = os.path.join('output', "bpm-{0}-{1}.fits".format(today_date, today_time))
    astropy.io.fits.writeto(filename=output_filename, data=combined_image_mask, header=astropy.io.fits.Header(header_dict))

    finished_time = time.time()

    time_diff = round(finished_time - current_time, 3)
    my_logger.info("Time elapsed: {0}".format(time_diff))


def get_image_hdus(absolute_path_to_images, suffix_list, hdu_index=0,):
    """
    Takes an absolute directory that contains images and converts those images to a list of header data objects.

    :param absolute_path_to_images: The absolute path of where the fits files are located. The FITS files must\
    be unpacked for this to work.

    :param hdu_to_retrieve: the index of the HDU you'd like. This will need to change to handle multi-extension fits files

    :param suffix_list: A list of file suffixes to use when retrieving the fits files. The suffixes should match\
    the filename format

    :return: A list of HDU objects

    """
    # Uses the suffix_list to construct a regex/wildcard pattern to search for files. So if your suffix list is:
    # a, b, c, this expression will be something like: '*[a|b|c].fits'
    if len(suffix_list) == 0:
        raise ValueError("No valid list of suffixes were passed.")

    file_wildcard_pattern = os.path.join(absolute_path_to_images, "*[{0}].fits".format('|'.join(suffix_list)))

    fits_images_list = [os.path.abspath(path) for path in glob.glob(file_wildcard_pattern)]

    hdu_list_for_images = [(astropy.io.fits.open(fits_image, mode='readonly'))[hdu_index] for fits_image in fits_images_list]

    return hdu_list_for_images


def sort_images_by_type(image_type, hdu_list):
    """
    Takes in a list of HDU objects, and only preserves the HDU objects whose obstype match the specified image type.

    :param image_type: Typically will be one of either, "BIAS", "DARK", or "SKYFLAT"
    :param hdu_list: A list of HDU objects
    :return: A list of HDU objects
    """

    # Should this function be condensed into a one-line list comprehension?
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
    global my_logger
    my_logger = logger.setup_custom_logger()
    my_logger.setLevel(logging.DEBUG)

    my_logger.debug("Starting main script.")

    if sys.argv[1] is None:
        raise TypeError("Missing parameter. No directory for source images was specified.")

    main(sys.argv[1])

