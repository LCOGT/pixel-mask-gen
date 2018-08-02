import fractions
import collections
import numpy
import astropy.stats
import logging


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


def test_adjacent_pixels(bad_pixel_list):
    """Test if adjacent pixels were marked as 'bad', this indicates some irregular activity, since the probability of \
    this happening naturally is very low.

    :param indiv_pixel_mask: A list of coordinates where each coordinate was the pixel that was marked as bad
    :return: The number of bad pixels that were adjacent to each other

    """
    max_neighboring_bad_pixels = 0

    for (row, col) in bad_pixel_list:
        # count the number of bad pixels that are adjacent to each other -- excluding diagonal
        adjacent_bad_pixel_count = sum((r, c) in bad_pixel_list for (r,c) in [(row,col-1), (row,col+1), (row-1, col),
                                                                              (row+1,col)])
        if adjacent_bad_pixel_count > max_neighboring_bad_pixels:
            max_neighboring_bad_pixels = adjacent_bad_pixel_count

    return max_neighboring_bad_pixels


def biases_processing(image_objects, sigma_min=7, sigma_max=7, pct_threshold=0.30):
    """**Algorithm**

    Compute the center quarter of the image, and then compute the median :math:`m` of the center quarter.

    Subtract every pixel in the median by :math:`m`.\

    Then perform the sigma clipping algorithm on the resulting pixels, and flag any that do not pass. Once you've\
    collected the list of pixels that don't pass, ensure that the same pixel appears at least 30% of the time before \
    marking it as truly 'bad'.

    :param image_objects: A list of image objects
    :param sigma_min: the lower sigma threshold to use
    :param sigma_max: the upper sigma threshold to use
    :return: A list of tuples where each tuple contains a pixel location that was flagged from the bias images.
    :rtype: list
    """
    corrected_images_list = []

    images_datas = [image.get_image_data() for image in image_objects]
    masked_indices_list = []
    filtered_std_deviation_list = []

    for image_data in images_datas:
        center_quarter = extract_center_fraction_region(image_data, fractions.Fraction(1, 4))
        center_quarter_median = numpy.median(center_quarter)

        corrected_image = numpy.subtract(image_data, center_quarter_median)

        corrected_images_list.append(corrected_image)

        sclipped_image = astropy.stats.sigma_clip(data=corrected_image, sigma_lower=sigma_min, sigma_upper=sigma_max, iters=5)

        masked_indices = numpy.transpose(numpy.ma.getmask(sclipped_image).nonzero())

        filtered_image = numpy.ma.masked_array(corrected_image, mask=numpy.logical_not(sclipped_image))

        masked_indices_list.append(masked_indices)

        stddev = numpy.std(filtered_image)

        filtered_std_deviation_list.append(stddev)

    # Return a flattened list containing all coordinates that failed the sigma clipping
    combined_list_of_failed_pixels = [tuple(coords) for sublist in masked_indices_list for coords in sublist]

    # the array containing the median-subtracted values that passed the image filtering. the sclipped image contains a
    # TRUE in every pixel that was removed, but since we want the pixels that were NOT removed, we invert the mask

    # once  you have the flattened list, count the frequencies of each pixel (i.e., how many times that specific pixel
    # appears in the list. Only use pixels who appear more than 30% of the time
    bias_bad_pixel_counter = collections.Counter(combined_list_of_failed_pixels)

    thresholded_bias_bad_pixel_list = [key for (key, value) in bias_bad_pixel_counter.items() if \
                                       (value >= (pct_threshold * len(images_datas)))]

    return thresholded_bias_bad_pixel_list


def darks_processing(image_objects, sigma_threshold=7, pct_threshold=0.10):
    r"""**Algorithm**

    1. Locate the overscan region of each image, and divide the entire image by the median of the overscan region.

    2. Then, divide the exposure time of each image, to obtain the dark current, in units of electrons/second.

    3. Filter any images that have a dark current of more than 35 electrons/second.

    :param image_objects: A list of image objects. See `image_objects.py` file.

    :return:  A list of tuples were each tuple contains a pixel location that was flagged from the darks images.

    :rtype: list
    """

    masked_indices_list = []

    logging.info("Beginning darks processing with {0} images".format(len(image_objects)))

    pixels_to_consider = [(253, 615), (363, 159), (733, 801), (1836, 2389), (1109, 2900), (955, 372), (1492, 2036), (1817, 838), (1627, 1164), (435, 216), (975, 2310), (2024, 818), (475, 1861), (544, 2844), (1782, 574), (1169, 1043), (763, 2055), (616, 2279), (692, 1615), (1250, 2094), (434, 1615), (38, 220), (1981, 2837), (278, 2254), (1834, 1176), (1642, 1880), (14, 1263), (1832, 2994), (1500, 1369), (1501, 48)]

    for image in image_objects:
        # Divide every pixel in the image by its exposure time, then store the new 'image' in a list

        bias_sect = image.get_image_header(key='BIASSEC')
        image_data = image.get_image_data()
        overscan_region_median = numpy.median(image_data[1:2048, 3100:3135])

        image_data = numpy.subtract(image_data, overscan_region_median)

        exposure_time = image.get_image_header(key='EXPTIME')

        image_data /= exposure_time

        filtered_image = numpy.ma.masked_less(image_data, 35)

        masked_indices = numpy.transpose(filtered_image.nonzero())

        masked_indices_list.append(masked_indices)

        for pixel in pixels_to_consider:
            print("pixel: {0}, electrons/second: {1}".format(pixel, image_data[pixel]))



    combined_list_of_bad_pixels =  [tuple(coords) for sublist in masked_indices_list for coords in sublist]

    darks_bad_pixel_counter = collections.Counter(combined_list_of_bad_pixels)

    thresholded_darks_bad_pixel_list = [key for (key, value) in darks_bad_pixel_counter.items() if \
                            (value >= (pct_threshold * len(image_objects)))]

    return thresholded_darks_bad_pixel_list


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

    logging.info("Beginning flats processing on {0} images".format(len(image_objects)))

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

    logging.info("Beginning center fraction extraction of image with shape: {0}".format(original_image_data.shape))
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
