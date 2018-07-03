import fractions
import numpy
import astropy.stats
import pdb
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

    mfiltered_array = astropy.stats.sigma_clip(data=image_array,\
                                    sigma_lower=sigma_low,\
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

def darks_processing(image_objects):
    r"""**Algorithm**

    Divide all dark images :math:`d_i` by their respective exposure time (:math:`e_i`) to get :math:`\bar{d_i}`

    Take the average of each image on a per-pixel basis, to get :math:`\bar{p}_{i,j}`

    Set some percentage threshold :math:`t` and remove all pixels outside the range :math:`\bar{p}_{i,j}  \pm t`

    :param image_objects:
    :return:  A list of tuples were each tuple contains a pixel location that was flagged from the darks images.
    :rtype: list
    """

    logging.info("Beginning darks processing with {0} images".format(len(image_objects)))
    corrected_image_list = []

    for image in image_objects:
        # Divide every pixel in the image by its exposure time, then store the new 'image' in a list
        exposure_time = image.get_image_header(key='EXPTIME')
        corrected_image = numpy.divide(image.get_image_data(), exposure_time)
        corrected_image_list.append(corrected_image)

    # A 3D array, where the x-y plane is the image plane and the z axis is the number of images
    total_image_cube = numpy.dstack(tuple(corrected_image_list))

    # compute the per-pixel mean, this should be in a 2D array
    per_pixel_mean = numpy.mean(total_image_cube, axis=2)

    # (scalar) median from every pixel
    pixels_median = numpy.median(per_pixel_mean)

    # any pixel that is outside of the range (median - 0.25 * median, median + 0.25 * median) gets masked
    scaling_factor = 0.25
    range_start, range_end = pixels_median - 0.25 * pixels_median, pixels_median + 0.25 * pixels_median

    filtered_array = numpy.ma.masked_outside(per_pixel_mean, range_start, range_end)

    masked_indices = numpy.transpose(numpy.ma.getmask(filtered_array).nonzero())

    return [tuple(coordinates) for coordinates in masked_indices.tolist()]


def flats_processing(image_objects):
    """**Algorithm**

    Compute the center quarter of the image, and then compute the median :math:`m` of the center quarter.

    Divide  every pixel in the image by :math:`m`.

    Store the data for each corresponding pixel in an array :math:`A_{i,j}`, and compute the standard deviation of the array,\
    for each value of :math:`i` and :math:`j`, and store this as :math:`\sigma_{A_{i,j}}`

    Take the median absolute deviation of all :math:`\sigma_{A_{i,j}}`, and flag any pixel location whose value is not within 3 standard \
    deviations of the median.


    :param image_objects: an array of image ojbects
    :return: A list of tuples where each tuple contains a pixel location that was flagged from the flats images.
    :rtype: list
    """

    logging.info("Beginning flats processing on {0} images".format(len(image_objects)))

    images_datas = [image.get_image_data() for image in image_objects]
    corrected_images_list = []

    for image_data in images_datas:
        center_quarter = extract_center_fraction_region(image_data, fractions.Fraction(1,4))
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

    # once you have the MAD, mask any values outside the range of 2*MAD?
    range_start, range_end = mad - 2*mad, mad + 2*mad

    filtered_array = numpy.ma.masked_outside(std_deviations_array, range_start, range_end)

    masked_indices = numpy.transpose(numpy.ma.getmask(filtered_array).nonzero())

    return [tuple(coordinates) for coordinates in masked_indices.tolist()]

def biases_processing(image_objects, sigma_min=5, sigma_max=5):
    """**Algorithm**

    Compute the center quarter of the image, and then compute the median :math:`m` of the center quarter.

    Subtract every pixel in the median by :math:`m`.\

    Then perform the sigma clipping algorithm on the resulting pixels, and flag any that do not pass.

    :param image_objects: A list of image objects
    :param sigma_min: the lower sigma threshold to use
    :param sigma_max: the upper sigma threshold to use
    :return: A list of tuples where each tuple contains a pixel location that was flagged from the bias images.
    :rtype: list
    """

    logging.info("Beginning processing on {0} bias images".format(len(image_objects)))

    corrected_images_list = []

    images_datas = [image.get_image_data() for image in image_objects]
    masked_indices_list = []

    for image_data in images_datas:
        center_quarter = extract_center_fraction_region(image_data, fractions.Fraction(1,4))
        center_quarter_median = numpy.median(center_quarter)

        corrected_image = numpy.subtract(image_data, center_quarter_median)

        corrected_images_list.append(corrected_image)

        sclipped_image = astropy.stats.sigma_clip(data=corrected_image,sigma_lower=sigma_min, sigma_upper=sigma_max)

        masked_indices = numpy.transpose(numpy.ma.getmask(sclipped_image).nonzero())

        masked_indices_list.append(masked_indices)

    # Return a flattened list containing all coordinates that failed the sigma clipping
    return [tuple(coords) for sublist in masked_indices_list for coords in sublist]


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


