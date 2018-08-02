import fractions
import collections
import numpy
import astropy.stats

import setup.LOGGER as logger


def sigma_clip_individual(image_array, sigma_hi, sigma_low):
    """Applies the median filter to one image, and returns back diagnostic information.

    :param image_array: A numpy array that represents the FITS file image
    :param sigma_hi: The upper limit in standard deviations to filter on
    :param sigma_low: The lower limit in standard deviations to filter on
    :returns mfiltered_array: a numpy array that has been median filtered (the filtered entries are masked)
    :returns masked_indices: the coordinates (x,y) that were masked by the median filter
    :returns percentage_mask: the percentage of pixels that were masked by the median filter

    """

    mfiltered_array = astropy.stats.sigma_clip(data=image_array,
                                        sigma_lower=sigma_low,
                                        sigma_upper=sigma_hi)

    masked_indices = numpy.transpose(numpy.ma.getmask(mfiltered_array).nonzero())

    percentage_masked = (len(masked_indices) / mfiltered_array.size) * 100

    return (mfiltered_array, masked_indices, percentage_masked)


def apply_bias_processing(hdu_objects, sigma_min=7, sigma_max=7, pct_threshold=0.30):
    """**Algorithm**

    Compute the center quarter of the image, and then compute the median :math:`m` of the center quarter.

    Subtract every pixel in the median by :math:`m`.\

    Then perform the sigma clipping algorithm on the resulting pixels, and flag any that do not pass. Once you've\
    collected the list of pixels that don't pass, ensure that the same pixel appears at least 30% of the time before \
    marking it as truly 'bad'.

    :param hdu_objects:

    """

    logger.debug("Beginning bias processing on {0} images".format(len(hdu_objects)))

    images_datas = [hdu.data for hdu in hdu_objects]

    list_of_masks = []

    for image_data in images_datas:
        center_quarter = extract_center_fraction_region(image_data, fractions.Fraction(1, 4))
        center_quarter_median = numpy.median(center_quarter)

        corrected_image = numpy.subtract(image_data, center_quarter_median)

        sclipped_image = astropy.stats.sigma_clip(data=corrected_image, sigma_lower=sigma_min, sigma_upper=sigma_max, iters=5)

        sclipped_image_as_masked_array = numpy.ma.getmaskarray(sclipped_image)

        list_of_masks.append(sclipped_image_as_masked_array)

    logger.debug("Completed bias processing")
    return list_of_masks


def extract_coordinates_from_header_string(header_string):
    """
    Utility for converting a specified string in the header into integers that can be used to slice an array.

    Example:  '[3100:3135,1:2048]' --> [3100, 3135, 1, 2048]


    :param header_string:
    :return: coordinate_list
    """

    two_d_coordinates = header_string.split(',')
    col_start, col_end = two_d_coordinates[0].split(':')
    row_start, row_end = two_d_coordinates[1].split(':')

    return [row_start, row_end, col_start, col_end]

def darks_processing(hdu_objects,dark_current_threshold=35):
    r"""**Algorithm**

    1. Locate the overscan region of each image, and divide the entire image by the median of the overscan region.

    2. Then, divide the exposure time of each image, to obtain the dark current, in units of electrons/second.

    3. Filter any images that have a dark current of more than 35 electrons/second.


    """

    logger.debug("Beginning darks processing on {0} images".format(len(hdu_objects)))

    list_of_masks = []

    for hdu in hdu_objects:
        # Divide every pixel in the image by its exposure time, then store the new 'image' in a list

        bias_section_header_string = hdu.header['BIASSEC']

        image_data = hdu.data

        overscan_region_coordinates = extract_coordinates_from_header_string(bias_section_header_string)

        # Remember that numpy does indexing in reverse order, e.g. (row : col) -> (col, row)
        cropped_image_data = image_data[overscan_region_coordinates[2] : overscan_region_coordinates[3],
                                        overscan_region_coordinates[1] : overscan_region_coordinates[4]]

        cropped_image_data_median = numpy.median(cropped_image_data)

        image_data -= cropped_image_data_median

        exposure_time = float(hdu.header['EXPTIME'])

        image_data /= exposure_time

        masked_image = image_data < dark_current_threshold

        list_of_masks.append(masked_image)

    logger.debug("Finished darks processing.")
    return list_of_masks


def flats_processing(image_objects, sigma_threshold=7):
    r"""**Algorithm**

    Compute the center quarter of the image, and then compute the median :math:`m` of the center quarter.

    Divide  every pixel in the image by :math:`m`.

    Store the data for each corresponding pixel in an array :math:`A_{i,j}`, and compute the standard deviation of the array,\
    for each value of :math:`i` and :math:`j`, and store this as :math:`\sigma_{A_{i,j}}`

    Take the median absolute deviation of all :math:`\sigma_{A_{i,j}}`, and flag any pixel location whose value is not within\
    :math:`s` standard deviations of the median.

    Recall that

    .. math::
    \sigma = k \cdot MAD

    Where for a normal distribution, :math:`k\approx 1.4826`.

    :param image_objects: an array of image ojbects
    :return: A list of tuples where each tuple contains a pixel location that was flagged from the flats images.
    :rtype: list
    """

    logger.debug("Beginning flats processing on {0} images".format(len(image_objects)))

    images_datas = [image.get_image_data() for image in image_objects]
    corrected_images_list = []

    for image_data in images_datas:
        center_quarter = extract_center_fraction_region(image_data, fractions.Fraction(1, 4))
        center_quarter_median = numpy.median(center_quarter)
        corrected_image = numpy.divide(image_data, center_quarter_median)
        corrected_images_list.append(corrected_image)

    # Create a 3D array out of all the corrected images, where (x,y) plane is the original image and the z-axis is what
    # connects the images
    corrected_images_cube = numpy.dstack(tuple(corrected_images_list))

    # Take the standard deviation of each pixel across all images
    # remember that axes are 0-indexed
    std_deviations_array = numpy.std(corrected_images_cube, axis=2)

    mad = astropy.stats.median_absolute_deviation(std_deviations_array)

    # once you have the MAD, mask any values outside the range of the ssigma threshold
    k = 1.4826
    range_start, range_end = mad - ((k * mad ) * sigma_threshold),\
                             mad + ((k * mad) * sigma_threshold)

    filtered_array = numpy.ma.masked_outside(std_deviations_array, range_start, range_end)

    masked_indices = numpy.transpose(numpy.ma.getmask(filtered_array).nonzero())


    return [tuple(coordinates) for coordinates in masked_indices.tolist()]


def extract_center_fraction_region(original_image_data, fraction):
    r"""Extract and return the center fraction of an image

    :param fraction: A Fraction object
    :param original_image_data:
    :return: The center fraction of the original image
    :rtype: numpy.ndarray

    """

    logger.info("Beginning center fraction extraction of image with shape: {0}".format(original_image_data.shape))
    row,col = original_image_data.shape

    if (row == 0) or (col == 0):
        raise ValueError('The array to be reduced has an invalid size.')

    row_start, row_end = int(((1 - fraction) * row) / 2), int((((1-fraction) * row) / 2) + row/2)

    if row == col: # the image is a square
        col_start, col_end = row_start, row_end

    else:
        col_start, col_end = int(((1 - fraction) * col)) / 2, int((((1-fraction) * col) / 2) + col/2)

    new_image_x, new_image_y = original_image_data.shape[0] // fraction.denominator, original_image_data.shape[1] // fraction.denominator

    extracted_image = original_image_data[new_image_y : -new_image_y, new_image_x : -new_image_x]

    return extracted_image
