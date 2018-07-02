# External
import astropy.io.fits
import astropy.stats
import numpy
import collections
import os
import yaml
import time
import datetime
import logging
import sys
import errno
import pdb


# Internal
import image_object
import image_processing
import fits_utilities

def main(arg1):
    """This function will execute when the module is run in main context.
     Loops over a range of camera IDs (0-99, must be two digits) and then tries to perform bad pixel masking on the
     corresponding camera with the specified camera prefix.

    :param arg1: the location of the config file, equivalent to sys.argv[1]: $python [arg0 = filename] \
    [arg1 = configuration file]

    :return: None.
    """
    empty_folder_count = 0
    for identifier in range(0, 99):
        camera_id_number = str(identifier).zfill(2)
        try:
            '''
            image_list, prefixes_list, full_camera_name = retrieve_image_directory_information(arg1,
                                                                                               camera_id_number)

            bias_array, dark_array, flat_array, image_header, rows, columns  = extract_data_from_files(image_list,prefixes_list)

            clean_bias_array, clean_dark_array, clean_flat_array = run_sigma_clipping([bias_array, dark_array,
                                                                                         flat_array], arg1)

            final_bpm_list = combine_bad_pixel_locations([clean_bias_array, clean_dark_array, clean_flat_array])

            final_bpm_mask = generate_mask_from_bad_pixels(final_bpm_list, rows, columns)

            output_to_FITS(final_bpm_mask, image_header, "{}_bpm.fits".format(full_camera_name))

            logger.info("Exiting main function.")
            '''

            image_list, prefixes_list, full_camera_name = retrieve_image_directory_information(arg1, camera_id_number)
            bias_array,dark_array,flat_array,image_header,rows,columns = extract_data_from_files(image_list, prefixes_list)

            bias_bad_pixels = image_processing.biases_processing(bias_array)

            dark_bad_pixels = image_processing.darks_processing(dark_array)

            flat_bad_pixels = image_processing.flats_processing(flat_array)

            final_bpm_list = combine_bad_pixel_locations([bias_bad_pixels, dark_bad_pixels, flat_bad_pixels])

            final_bpm_mask = generate_mask_from_bad_pixels(final_bpm_list, rows, columns)

            fits_utilities.output_to_FITS(final_bpm_mask, {}, "{}_bpm.fits".format(full_camera_name))

            logger.info('Exiting main function.')

        except FileNotFoundError:
            logger.info("Unable to find any folders that contained the desired camera prefix and identifier ({0})".format(camera_id_number))
            empty_folder_count += 1
            continue

    if empty_folder_count == 99:
        logger.warning("No folders matching ANY of the camera prefix and identifiers were found. Check this folder.")

def setup_custom_logger(name='pixel-mask-gen'):
    """This function defines the custom logging instance for use in the module. Logs will be outputted in the form:
    <YYYY-MM-DD HH:MM:SS> <LEVEL> <MSG>

    **Log Level Information**

    INFO : General information, replacement to print().

    WARNING : Critical warnings that require manual review.

    ERROR : Unrecoverable failures that cause exceptions.

    :param name: The name to use to define a logging instance. Defaults to the name of the repository: 'pixel-mask-gen'
    :return: A logger object that will be used throughout the program.
    :rtype: logging.Logger

    """
    formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)
    handler.setLevel(logging.DEBUG)
    logger = logging.getLogger(name)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    return logger

def retrieve_image_directory_information(config_file_location, camera_identifier):
    """Parses the YAML configuration file to retrieve information about where the images will be located

    :param config_file_location: The relative path of the YAML configuration file

    :returns image_list: A list of absolute file paths representing the FITS images to filter.\

    :returns prefixes_list: A list of prefixes corresponding to image calibration type (bias, dark, light).

    :return: the camera designation prefix (letters) plus id code. Will be used to search for source image folders and \
    will be used to name the output fits file (since one bad pixel mask matches up with every camera)

    """

    if (os.path.exists(config_file_location)):
        logger.info('Opening configuration file located at {}'.format(config_file_location))

    else:
        logger.error('Unable to locate file at {}'.format(config_file_location))
        #logger.error(msg)
        #raise OSError(msg)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), config_file_location)

    with open(config_file_location) as yaml_file:
        # Store as dictionary
        configData = yaml.safe_load(yaml_file)

        if (configData is None):
            raise ValueError('Yaml file exists but appears to be empty.')

        if (len(configData) <= 1):
            raise ValueError('Yaml file exists but was incorrectly parsed, please check that it is valid.')

        directory_info = configData['directories']

        # Get the latest date string to use.
        date_string = parse_config_file_date(directory_info)

        # Check that top level directory is a real directory before you even start. If its not, then exit -- if the top directory isn't valid then
        # none of the subdirectories can be, so don't waste your time checking them...
        if os.path.isdir(directory_info['top_directory']):
            logger.info("Top directory ({0}) appears valid, will continue searching for images.".format(directory_info['top_directory']))

        else:
            #raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), directory_info['top_directory'])
            sys.exit(str(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), directory_info['top_directory'])))

        # Locate all the image files and how many there are
        image_folder = os.path.join(directory_info['top_directory'], directory_info['camera_prefix'] + camera_identifier, date_string)
        # Once you have the image folder, check if it exists

        logger.info("Searching for images in the path: {0}".format(image_folder))

        if os.path.isdir(image_folder):
            logger.info("Image folder is {}".format(image_folder))

        else:
            #msg = "The image folder you requested: {} does not appear to exist.".format(image_folder)
            #raise OSError("The image folder you requested does not appear to exist -- will search for another camera prefix.")
            logger.error("Image folder you requested does not appear to exist, throwing FileNotFoundError")
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), image_folder)

        # Stores the absolute path of the images in a list
        image_list = [os.path.abspath(os.path.join(image_folder, filename)) for filename in os.listdir(image_folder) if \
                      os.path.isfile(os.path.abspath(os.path.join(image_folder, filename)))]
        prefixes_list = [directory_info['bias_prefix'], directory_info['dark_prefix'], directory_info['flats_prefix']]

        if len(image_list) < 1:
            msg = "The image folder exists but appears to be empty."
            logger.error(msg)
            raise ValueError("The image folder exists but appears to be empty.")

        else:
            logger.info("{} images found".format(len(image_list)))

        return (image_list, prefixes_list, directory_info['camera_prefix'] + camera_identifier)


def extract_data_from_files(image_list, prefixes_array):
    """Takes in a list of FITS filenames, converts FITS data into numpy arrays, and then stores that FITS info according \
    to what calibration type the image corresponds to.

    :param image_list: A list of absolute filenames for images to search for.

    :param prefixes_array: An array of prefixes from the config.yml file that denote bias,dark,light

    :returns (bias|dark|flat) array: a list of numpy arrays, each numpy array represents one image \
    rows and cols are the sizes of the images

    :returns image_header: the image header stored from the last image

    """

    f_array, d_array, b_array = [], [], []

    logger.info("Image list received: {0}".format(image_list))
    logger.info("Corresponding prefix list received: {0}".format(prefixes_array))

    if len(image_list) < 1:
        msg = "No images were passed in, unable to extract data."
        #logger.error(msg)
        raise ValueError(msg)

    for image_filename in image_list:
        # check which prefix the filename ends with, then, the first letter of the matching prefix
        # will tell you what array it should go into
        for prefix in prefixes_array:
            if image_filename.endswith("{0}.fits".format(prefix)):
                image_data, image_headers, image_shape = fits_utilities.read_individual_fits_file(image_filename)

                # Create image object
                image = image_object.ImageObject(image_data, image_headers)

                arr_string = prefix[0] + "_array"
                #eval(arr_string + '.append(image_data)') # example: b_array.append(image_data)
                eval(arr_string + '.append(image)') # example: b_array.append(image)

                # once you've found the matching prefix, stop checking for other prefixes -- at most once prefix can match
                break

    # not all arrays can be empty here or something went wrong
    if ( len(b_array) == len(d_array) == len(f_array) == 0 ):
        raise ValueError('All image calibration arrays are empty.')

    else:
        logger.info("Data extraction from FITS files complete. {0} images parsed in total".format(len(b_array) + \
                                                                                                  len(f_array) + \
                                                                                                  len(d_array)))
    # return the arrays in alphabetical order
    return (b_array, d_array, f_array, image_headers, image_shape[0], image_shape[1])

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


def run_sigma_clipping(images_arrays, config_file_location):
    """Once you have a list of all the images (where each image is a numpy array), you are now ready to
    apply the proper median filtering on them.

    :param images_arrays: A list of three lists (bias, dark, flat) where each array contains a list of numpy arrays \
    that each represent one image
    :return: A list of lists, where each sublist contains all pixels that were marked as 'known-bad'

    """
    # Run median filtering on each type of image: bias, dark, flats. Store the arrays of which indices
    # were masked. At the end, use the indices that appeared in EVERY mask (e.g. only mark a pixel
    # as bad if it was filtered in X percentage of images)

    with open(config_file_location) as yaml_file:
        configData = yaml.safe_load(yaml_file)
        statistics_info = configData['statistics']
        sigma_hi = statistics_info['sigma_high']
        sigma_low = statistics_info['sigma_low']
        threshold = statistics_info['threshold']
        max_pct_removed = statistics_info['max_pct_removed']

    # Each object below is list of arrays, where each array will contains the locations of bad pixels
    dark_bpm_collection, flat_bpm_collection, bias_bpm_collection = [], [], []

    for index_cal_array, cal_array in enumerate(images_arrays):
        for index_image_array, image_array in enumerate(cal_array):

            filtered_array, masked_indices, percentage_masked = image_processing.sigma_clip_individual(image_array, sigma_hi, sigma_low)

            logger.info("{0} ({1}% of total) pixels failed the sigma clipping on image #{2}".format(len(masked_indices),
                                                                                                    percentage_masked,
                                                                                                    index_image_array))

            logger.info("Failing pixels occured at: {0}".format(masked_indices))

            # make sure too many pixels arent getting filtered, but if they are, increase the minimum std dev that you
            # allow
            while (percentage_masked >= max_pct_removed):
                logger.info("Too many pixels failed sigma clipping. Increasing sigma values from {0} to {1} (low) and \
                from {2} to {3} (high)".format(sigma_low, sigma_low +1, sigma_hi, sigma_hi + 1))
                sigma_hi += 1
                sigma_low += 1
                filtered_array, masked_indices, percentage_masked = image_processing.sigma_clip_individual(image_array, sigma_hi,sigma_low)
                logger.info("{0} ({1}% of total) pixels failed the sigma clipping on image #{2}".format(len(masked_indices),
                                                                                                        percentage_masked,
                                                                                                        index_image_array))
            if (index_cal_array == 0):
                bias_bpm_collection.append(masked_indices)

            elif (index_cal_array == 1):
                dark_bpm_collection.append(masked_indices)

            elif (index_cal_array == 2):
                flat_bpm_collection.append(masked_indices)

    # once you have all the bad pixel indices from each image, separated by calibration type,see how often each bad
    # pixel appears, but to count frequency of appearance across images, you'll need to flatten the list
    # a 1D array storing every coordinate that is marked as 'bad', in tuple form
    bias_bpm_flattened, dark_bpm_flattened, flat_bpm_flattened = [], [], []
    '''
    # TODO: Remove code reuse
    for sublist in bias_bpm_collection:
        for coords in sublist:
            bias_bpm_flattened.append(tuple(coords))

    for sublist in dark_bpm_collection:
        for coords in sublist:
            dark_bpm_flattened.append(tuple(coords))

    for sublist in flat_bpm_collection:
        for coords in sublist:
            flat_bpm_flattened.append(tuple(coords))
    '''

    bias_bpm_flattened = generate_flattened_list(bias_bpm_collection)
    dark_bpm_flattened = generate_flattened_list(dark_bpm_collection)
    flat_bpm_flattened = generate_flattened_list(flat_bpm_flattened)

    bias_bpm_counter = collections.Counter(bias_bpm_flattened)
    dark_bpm_counter = collections.Counter(dark_bpm_flattened)
    flat_bpm_counter = collections.Counter(flat_bpm_flattened)


    # The threshold amount is threshold_pct * length of each array, any pixels appearing with threshold equal to or
    # higher than that are marked as bad
    # TODO: Remove the code reuse here
    thresholded_bias_array = [key for (key, value) in bias_bpm_counter.items() if (value >= (float(threshold)/100 * len(bias_bpm_collection)))]
    thresholded_dark_array = [key for (key, value) in dark_bpm_counter.items() if (value >= (float(threshold)/100 * len(dark_bpm_collection)))]
    thresholded_flat_array = [key for (key, value) in flat_bpm_counter.items() if (value >= (float(threshold)/100 * len(flat_bpm_collection)))]

    return (thresholded_bias_array, thresholded_dark_array, thresholded_flat_array)


def generate_mask_from_bad_pixels(bad_pixel_location_array, x_dimension, y_dimension):
    """The array in the parameter given has the coordinates of all known-bad pixels. Use these to create a mask \
    where the masks have True for bad pixels and false otherwise

    :param bad_pixel_location_array: An array containing the coordinates (stored as tuples) of all the known-bad pixels

    :return masked_array: a representation of the bad pixel coordinates a boolean array, where known-bad pixels are\
    marked by true and safe pixels are marked as false

    """

    masked_array = numpy.zeros((x_dimension, y_dimension), dtype=bool)
    size = x_dimension * y_dimension
    logger.info("Final count: {0} ({1}% of total) bad pixels detected.".format(len(bad_pixel_location_array),
                                                                               (len(bad_pixel_location_array)/size)*100))
    for coordinates in bad_pixel_location_array:
        masked_array[coordinates] = True

    return masked_array

def combine_bad_pixel_locations(arrays_of_bad_pixels):
    """Given an array of 3 arrays, (where each inner array contains a list of pixels that are known bad), take the set union

    :param arrays_of_bad_pixels: 3 arrays that contain a list of pixels that were set as known-bad during the filtering \
    process

    :return: a *set* that describes every bad pixel that appeared, across all calibration types --\
    including bias, dark, and flats

    """

    for index, arr in enumerate(arrays_of_bad_pixels):
        logger.info("{0} bad pixels were detected for calibration type #{1}".format(len(arr), index))
        indiv_file_path = os.path.join("debug","{0}_bpm.txt".format(index))
        logger.info('Writing bad pixel to debugging directory')
        with open(indiv_file_path, 'w') as output:
            output.write(''.join(map(str, arr)))
            output.close()

        # If the array contains bad pixels that are adjacent, flag it
        neighboring_bad_pixels = image_processing.test_adjacent_pixels(arr)
        msg = "{0} neighboring bad pixels were found.".format(neighboring_bad_pixels)
        if neighboring_bad_pixels < 1:
            logger.info(msg)

        elif neighboring_bad_pixels == 1:
            logger.warning(msg + "Potentially review images.")

        else:
            logger.error(msg + "Review images now, image mask is corrupt")

    final_bad_pixel_set = set()

    for subarray in arrays_of_bad_pixels:
        #pdb.set_trace()
        for coords in subarray:
            final_bad_pixel_set.add(tuple(coords))

    # As a sanity check, print the list of coordinates into a text file for later examination
    final_file_path = 'debug/combined_bpm_list.txt'
    with open(final_file_path, 'w') as output:
        output.write(''.join(map(str, final_bad_pixel_set)))
        output.close()

    return final_bad_pixel_set

def get_last_date_of_unit(todays_date, string):
    """Returns the latest date that falls within the last month or week.
    For example, if today is Thursday, June 7th, 2018, then get_last_date_of_unit(week) will return \
    the end of last week, which was Saturday, June 2nd, 2018.

    :param todays_date: A datetime object representing today's date.
    :param string: 'month','year','week'
    :return: A datetime string

    """
    if string == 'week':
        return str(todays_date - datetime.timedelta(days=(todays_date.isoweekday() + 1)))

    elif string == 'month':
        return str(todays_date - datetime.timedelta(days=(todays_date.day + 1)))

    elif string == 'year':
        return str(datetime.date(year=(todays_date.year - 1), month=12, day=31))

def parse_config_file_date(directory_info):
    """Takes in a dictionary object representing the information from the YML file for configuration and returns\
    a string representation of the date that was requested.

    :param directory_info: A dictionary representation of the 'directories' node in the YML file
    :return: An ISO-8601 string representation (YYYMMDD) that will tell the program what directory to use

    """
    VALID_RELATIVE_DATES = ['yesterday', 'last_month', 'last_week']

    logger.info("Directory info received was: {0}".format(directory_info))

    if directory_info['date']['exact'] != False:
        # Convert the datetime string to a datetime object, and see if theres an error
        try:
            time.strptime(directory_info['date']['exact'], '%Y%m%d')
        except ValueError:
            logger.error("Error parsing the exact date given from configuration file. Reverting to yesterday's date.")

        finally: # Revert to yesterday's date
            return (datetime.datetime.today() - datetime.timedelta(days=1)).strftime('%Y%m%d')

    else:
        if directory_info['date']['relative'] in VALID_RELATIVE_DATES:
            my_date = directory_info['date']['relative']
            if my_date == 'yesterday':
                return (datetime.datetime.today() - datetime.timedelta(days=1)).strftime('%Y%m%d')

            elif (my_date == 'last_week'): # Gets the latest date that falls last week
                return get_last_date_of_unit(datetime.datetime.today().strftime('%Y%m%d'), 'week')

            elif (my_date == 'last_month'):
                return get_last_date_of_unit(datetime.datetime.today().strftime('%Y%m%d'), 'month')

        else:
            logger.info("Invalid relative date in configuration file, reverting to yesterday's date.")
            return (datetime.datetime.today() - datetime.timedelta(days=1)).strftime('%Y%m%d')

global logger
logger = setup_custom_logger('pixel-mask-gen')

if __name__ == '__main__':
    main(sys.argv[1])

