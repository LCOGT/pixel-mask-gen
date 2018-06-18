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
import functools
import pdb

def setup_custom_logger(name):
    formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)
    handler.setLevel(logging.DEBUG)
    logger = logging.getLogger(name)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    #logger.addHandler(handler)

    return logger
def locate_image_list():
    """

    :return:
    """
    with open('config.yml') as yaml_file:
        # Store as dictionary
        configData = yaml.safe_load(yaml_file)
        directory_info = configData['directories']

        # Get the latest date string to use.
        date_string = parse_config_file_date(directory_info)

        # Locate all the image files and how many there are.
        # the format is camera/date/<start> or camera/date/<end>
        image_folder = os.path.join(directory_info['top_directory'], directory_info['camera_prefix'],
                                    date_string)
        #print("image_folder", image_folder)
        logger.info("Image folder is {}".format(image_folder))

        # Stores the absolute path of the images in a list
        #image_list = [os.path.abspath(file) for file in os.listdir(image_folder)]
        image_list = [os.path.abspath(os.path.join(image_folder, filename)) for filename in os.listdir(image_folder)]
        #print("image list:")

        if len(image_list) < 1:
            raise ValueError('No images found -- please check that path was correct.')

        else:
            logger.info("{} images found".format(len(image_list)))

        return image_list

def extract_data_from_files(image_list, prefixes_array):
    """
    Takes in a list of FITS file images, returns 3 lists of lists of numpy arrays (bias, flat, dark)
    where each array contains the data stored in the image (and returns the header for each file)
    :param image_list:
    :param prefixes_array: An array of prefixes from the config.yml file that denote bias,dark,light
    :return:
    """

    f_array, d_array, b_array = [], [], []

    for image_filename in image_list:
        # check which prefix the filename ends with, then, the first letter of the matching prefix
        # will tell you what array it should go into
        for prefix in prefixes_array:
            if image_filename.endswith("{0}.fits".format(prefix)):
                image_file = astropy.io.fits.open(image_filename)
                image_data = image_file[0].data
                rows, cols = image_data.shape
                if (rows == 0) or (cols == 0):
                    logger.error("The image {0} has no data".format(image_filename))
                    raise ValueError("The image you are trying to read has no data!")

                else:
                    logger.info("Located image: {0} having shape: {1} by {2}".format(image_filename, rows, cols))
                    image_file.close()
                    #print("x_size:", x_size, "y_size:", y_size, "for image: ", image_filename)
                    arr_string = prefix[0] + "_array"
                    eval(arr_string + '.append(image_data)') # example: b_array.append(image_data)

                    # once you've found the matching prefix, stop checking for other prefixes -- at
                    # most once prefix can match
                    break

    return (b_array, d_array, f_array) # return this in alphabetical order

def filter_individual(image_array, sigma_hi, sigma_low):
    mfiltered_array = astropy.stats.sigma_clip(data=image_array,\
                                    sigma_lower=sigma_low,\
                                    sigma_upper=sigma_hi)

    masked_indices = numpy.transpose(numpy.ma.getmask(mfiltered_array).nonzero())

    percentage_masked = (len(masked_indices) / mfiltered_array.size) * 100

    return (mfiltered_array, masked_indices, percentage_masked)


def run_median_filtering(images_arrays):
    """
    Once you have a list of all the images (where each image is a numpy array), you are now ready to
    apply the proper median filtering on them
    :param images_arrays: A list of three arrays (b, d, f) where each array contains a list of numpy arrays
    :return: A list of pixel indices (x,y) that contain bad pixels as defined by a certain threshold
    """

    """
    Run median filtering on each type of image: bias, dark, flats. Store the arrays of which indices 
    were masked. At the end, use the indices that appeared in EVERY mask (e.g. only mark a pixel 
    as bad if it appeared bad in X percentage of images)
    """
    with open('config.yml') as yaml_file:
        configData = yaml.safe_load(yaml_file)
        statistics_info = configData['statistics']
        sigma_hi = statistics_info['sigma_high']
        sigma_low = statistics_info['sigma_low']
        threshold = statistics_info['threshold']
        max_pct_removed = statistics_info['max_pct_removed']

    # Run sigma clipping on each image, see what points were masked

    # Each object below is list of arrays, where each array contains the locations of bad pixels
    dark_bpm_collection, flat_bpm_collection, bias_bpm_collection = [], [], []

    for index_cal_array, cal_array in enumerate(images_arrays):
        for index_image_array, image_array in enumerate(cal_array):
            """
            mfiltered_array = astropy.stats.sigma_clip(data=image_array,\
                                    sigma_lower=sigma_low,\
                                    sigma_upper=sigma_hi)
            """

            # figure out the indices of the pixels that were masked for each image, and add them to
            #  the corresponding bpm collection

            #masked_indices = numpy.transpose(numpy.ma.getmask(mfiltered_array).nonzero())

            filtered_array, masked_indices, percentage_masked = filter_individual(image_array, sigma_hi, sigma_low)

            #print(len(masked_indices), "pixels failed sigma clipping on image #", index_image_array, "of the current type")
            logger.info("{0} ({1}% of total) pixels failed the sigma clipping on image #{2}".format(len(masked_indices), percentage_masked, index_image_array))

            # make sure too many pixels arent getting filtered, but if they are, increase the minimum std dev that you
            # allow
            while (percentage_masked >= max_pct_removed):
                logger.info("Too many pixels failed sigma clipping. Lowering sigma values from {0} to {1} (low) and from {2} to {3} (high)".format(sigma_low, sigma_low +1, sigma_hi, sigma_hi + 1))
                sigma_hi += 1
                sigma_low += 1
                filtered_array, masked_indices, percentage_masked = filter_individual(image_array, sigma_hi,sigma_low)
                logger.info("{0} ({1}% of total) pixels failed the sigma clipping on image #{2}".format(len(masked_indices), percentage_masked, index_image_array))



            if (index_cal_array == 0):
                bias_bpm_collection.append(masked_indices)

            elif (index_cal_array == 1):
                dark_bpm_collection.append(masked_indices)

            elif (index_cal_array == 2):
                flat_bpm_collection.append(masked_indices)


    # once you have all the bad pixel indices from each image, separated by calibration type,
    # see how often each bad pixel appears, but to count frequency of appearance across images, you'll need to flatten
    # the list
    """
    bias_bpm_counter = collections.Counter(*bias_bpm_collection)
    dark_bpm_counter = collections.Counter(*dark_bpm_collection)
    flat_bpm_counter = collections.Counter(*flat_bpm_collection)
    """
    # a 1D array storing every coordinate that is marked as 'bad', in tuple form
    bias_bpm_flattened, dark_bpm_flattened, flat_bpm_flattened = [], [], []
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


    bias_bpm_counter = collections.Counter(bias_bpm_flattened)
    dark_bpm_counter = collections.Counter(dark_bpm_flattened)
    flat_bpm_counter = collections.Counter(flat_bpm_flattened)


    # The threshold amount is threshold_pct * length of each array, any pixels appearing with threshold
    # equal to or higher than that are marked as bad
    # TODO: Remove the code reuse here
    thresholded_bias_array = [key for (key, value) in bias_bpm_counter.items() if (value >= (float(threshold)/100 * len(bias_bpm_collection)))]
    thresholded_dark_array = [key for (key, value) in dark_bpm_counter.items() if (value >= (float(threshold)/100 * len(dark_bpm_collection)))]
    thresholded_flat_array = [key for (key, value) in flat_bpm_counter.items() if (value >= (float(threshold)/100 * len(flat_bpm_collection)))]

    #pdb.set_trace()

    return thresholded_bias_array, thresholded_dark_array, thresholded_flat_array


def generate_mask_from_bad_pixels(bad_pixel_location_array, x_dimension, y_dimension):
    """
    The array in the parameter given has the coordinates of all known-bad pixels. Use these to create a mask where the masks have True for bad pixels and false otherwise
    :param bad_pixel_location_arrays:
    :return:
    """
    """
    for index_bad_pixel_array, bad_pixel_array in enumerate(bad_pixel_location_arrays):
        masked_array = numpy.zeros((x_dimension, y_dimension), dtype=bool)
        #pdb.set_trace()
        #masked_array[bad_pixel_array] = True
        # Bad pixel array contains coordinates of all known bad pixels, so set these coordinates to true in the array
        #print(len(bad_pixel_array), "known-bad pixels found")
        logger.info("{0} ({1}% of total) bad pixels detected".format(len(bad_pixel_array), len(bad_pixel_array) / (x_dimension * y_dimension)))
        for coordinates in bad_pixel_array:
            masked_array[coordinates] = True

        array_masks.append(masked_array)

        if (debug == True):
            file_path = "debug/{0}_bpm.txt".format(index_bad_pixel_array)
            with open(file_path, 'a') as output:
                output.write(''.join(map(str, bad_pixel_array)))
            output.close()

    """

    masked_array = numpy.zeros((x_dimension, y_dimension), dtype=bool)
    size = x_dimension * y_dimension
    logger.info("Final count: {0} ({1}% of total) bad pixels detected.".format(len(bad_pixel_location_array), (len(bad_pixel_location_array) / size) * 100))
    for coordinates in bad_pixel_location_array:
        masked_array[coordinates] = True

    return masked_array


def combine_bad_pixel_locations(arrays_of_bad_pixels):
    """
    Given an array of 3 arrays, (where each inner array contains a list of pixels that are known bad), take the set union
    :param array_of_masks:
    :return:
    """
    for index, arr in enumerate(arrays_of_bad_pixels):
        logger.info("{0} bad pixels were detected for calibration type #{1}".format(len(arr), index))
        indiv_file_path = "debug/{0}_bpm.txt".format(index)
        with open(indiv_file_path, 'w') as output:
            output.write(''.join(map(str, arr)))
            output.close()
    #final_badpixellist_combined = functools.reduce(numpy.union1d, (arrays_of_bad_pixels[0], arrays_of_bad_pixels[1], arrays_of_bad_pixels[2]))

    final_bad_pixel_set = set()
    for subarray in arrays_of_bad_pixels:
        for coords in subarray:
            final_bad_pixel_set.add(coords)

    final_file_path = 'debug/combined_bpm.txt'
    with open(final_file_path, 'w') as output:
        output.write(''.join(map(str, final_bad_pixel_set)))
        output.close()

    return final_bad_pixel_set


def output_to_FITS(image_data, header_dict, filename):
    """
    Generates a FITS v4 file from image data.
    :param image_data: A numpy array, will be used for the primary header.
    :param header_dict: A dictionary used for the header data unit key/value pairs.
    :param filename: The destination file name for the FITS file.
    :return: None
    """
    # See: https://fits.gsfc.nasa.gov/fits_primer.html for info about FITS

    if image_data.size < 1:
        raise BaseException("No data to write to Primary HDU.")

    else:
        new_hdu = astropy.io.fits.PrimaryHDU(image_data.astype(int))
        #print("writing array of size", image_data.shape, "to fits file. Array contains", image_data.sum(), "bad pixels in mask")
        logging.info("Writing array of size {0} to .fits file; array contains {1} bad pixels in mask.".format(image_data.shape, image_data.sum()))
        new_hdu_list = astropy.io.fits.HDUList([new_hdu])


    for key, value in header_dict.items():
        new_hdu_list[0].header.set(key, value)

    new_hdu_list.writeto(filename, overwrite=True, output_verify='exception',checksum=True)

    new_hdu_list.close()


def get_last_date_of_unit(todays_date, string):
    """
    Returns the latest date that falls within the last month or week.
    For example, if today is Thursday, June 7th, 2018, then get_last_date_of_unit(week) will return
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
    """

    :param directory_info:
    :return:
    """
    VALID_RELATIVE_DATES = ['yesterday', 'last_month', 'last_week']

    if directory_info['date']['exact'] == 'true':
        # Convert the datetime string to a datetime object, and see if theres an error
        try:
            time.strptime(directory_info['date']['exact'], '%Y%m%d')
        except ValueError:
            # TODO: Replace with logging/warning
            print("Error parsing the exact date given from config.yml. Reverting to yesterday's date.")

        finally: # Revert to yesterday's date
            return (datetime.today() - datetime.timedelta(days=1)).strftime('%Y%m%d')

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
            # TODO: Replace with logging/warning
            logger.info("Invalid relative date in config.yml, reverting to yesterday's date.")
            return (datetime.today() - datetime.timedelta(days=1)).strftime('%Y%m%d')


if __name__ == '__main__':
    # run the sequence
    print("Sequence starting now.")
    logger = setup_custom_logger('pixel_mask_gen')

    image_list = locate_image_list()
    # manually hardcode prefix list for now?
    bias_array, dark_array, flat_array = extract_data_from_files(image_list, ['b00','f00','d00'])
    clean_bias_array, clean_dark_array, clean_flat_array = run_median_filtering([bias_array, dark_array, flat_array])

    # temporarily hardcode image size?
    #array_of_masks = generate_masks_from_bad_pixels([clean_bias_array, clean_dark_array, clean_flat_array], 2112, 3136)
    """
    for index, mask_array in enumerate(array_of_masks):
        output_to_FITS(mask_array, {}, str(index) + '_bpm.fits')
        
    """
    final_bpm_list = combine_bad_pixel_locations([clean_bias_array, clean_dark_array, clean_flat_array])
    final_bpm_mask = generate_mask_from_bad_pixels(final_bpm_list, 2112, 3136)
    output_to_FITS(final_bpm_mask, {}, "final_bpm.fits")


    print("Sequence ended.")


