import astropy.io.fits
import astropy.stats
import numpy
import collections
import os
import yaml
import time
import datetime
import pdb


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
        print("image_folder", image_folder)

        # Stores the absolute path of the images in a list
        #image_list = [os.path.abspath(file) for file in os.listdir(image_folder)]
        image_list = [os.path.abspath(os.path.join(image_folder, filename)) for filename in os.listdir(image_folder)]
        print("image list:")
        print(image_list)

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
                image_file.close()
                arr_string = prefix[0] + "_array"
                eval(arr_string + '.append(image_data)') # example: b_array.append(image_data)

                # once you've found the matching prefix, stop checking for other prefixes -- at
                # most once prefix can match
                break

    return (b_array, d_array, f_array) # return this in alphabetical order

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

    # Run sigma clipping on each image, see what points were masked

    # Each object below is list of arrays, where each array contains the locations of bad pixels
    dark_bpm_collection, flat_bpm_collection, bias_bpm_collection = [], [], []

    for index_cal_array, cal_array in enumerate(images_arrays):
        for image_array in cal_array:
            mfiltered_array = astropy.stats.sigma_clip(data=image_array,\
                                    sigma_lower=sigma_low,\
                                    sigma_upper=sigma_hi)

            # figure out the indices of the pixels that were masked for each image, and add them to
            #  the corresponding bpm collection

            masked_indices = numpy.transpose(numpy.ma.getmask(mfiltered_array).nonzero())

            if (index_cal_array == 0):
                bias_bpm_collection.extend(masked_indices)

            elif (index_cal_array == 1):
                dark_bpm_collection.extend(masked_indices)

            elif (index_cal_array == 2):
                flat_bpm_collection.extend(masked_indices)


    # once you have all the bad pixel indices from each image, separated by calibration type,
    # see how often each bad pixel appears
    pdb.set_trace()
    bias_bpm_counter = collections.Counter(*bias_bpm_collection)
    dark_bpm_counter = collections.Counter(*dark_bpm_collection)
    flat_bpm_counter = collections.Counter(*flat_bpm_collection)

    # The threshold amount is threshold_pct * length of each array, any pixels appearing with threshold
    # equal to or higher than that are marked as bad
    # TODO: Remove the code reuse here
    thresholded_bias_array = [value for (key, value) in bias_bpm_counter.items() if (value >= int(threshold) * len(bias_bpm_counter))]
    thresholded_dark_array = [value for (key, value) in dark_bpm_counter.items() if (value >= int(threshold) * len(dark_bpm_counter))]
    thresholded_flat_array = [value for (key, value) in flat_bpm_counter.items() if (value >= int(threshold) * len(flat_bpm_counter))]

    return thresholded_bias_array, thresholded_dark_array, thresholded_flat_array


def generate_mask_from_bad_pixels(bad_pixel_location_arrays, x_dimension, y_dimension):
    """
    Each array in the parameter given has the coordinates of all known-bad pixels. Use these to create 3 different masks (one for
    each calibration type) where the masks have True for bad pixels and false otherwise
    :param bad_pixel_location_arrays:
    :return:
    """
    array_masks = []

    for bad_pixel_array in bad_pixel_location_arrays:
        masked_array = numpy.zeros((y_dimension, x_dimension), dtype=bool)
        masked_array[bad_pixel_array] = True
        array_masks.append(masked_array)


    return array_masks


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
        new_hdu = astropy.io.fits.PrimaryHDU(image_data)
        new_hdu_list = astropy.io.fits.HDUList([new_hdu])

    for key, value in header_dict.items():
        new_hdu_list[0].header.set(key, value)

        new_hdu_list.writeto(filename, clobber=True, output_verify='exception')

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
            return str(datetime.today() - datetime.timedelta(days=1))

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
            print("Invalid relative date. Reverting to yesterday's date.")
            pass





if __name__ == '__main__':
    # run the sequence
    print("Sequence starting now.")

    image_list = locate_image_list()
    # manually hardcode prefix list for now?
    bias_array, dark_array, flat_array = extract_data_from_files(image_list, ['b00','f00','d00'])
    clean_bias_array, clean_dark_array, clean_flat_array = run_median_filtering([bias_array, dark_array, flat_array])

    # temporary hardcode x and y axis limits?
    array_of_masks = generate_mask_from_bad_pixels([clean_bias_array, clean_dark_array, clean_flat_array], 2000, 2000)
    for index, mask_array in array_of_masks:
        output_to_FITS(mask_array, {}, str(index) + '_bpm.fits')


