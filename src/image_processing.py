import numpy
import astropy.stats
import logging

my_logger = logging.getLogger(__name__)

# see: https://en.wikipedia.org/wiki/Median_absolute_deviation
global mad_constant
mad_constant = 1.4826


def apply_bias_processing(hdu_objects,
                          sigma_min=7,
                          sigma_max=7,
                          sigma_clip_iters=None,
                          sigma_threshold=7):
    r"""**Algorithm**

    Compute the center quarter of the image, and then compute the median :math:`m` of the center quarter.

    Subtract every pixel in the median by :math:`m`.\

    Then perform the sigma clipping algorithm on the resulting pixels, and flag any that do not pass. Once you've\
    collected the list of pixels that don't pass, ensure that the same pixel appears at least 30% of the time before \
    marking it as truly 'bad'.

    :param hdu_objects:

    """

    my_logger.debug("Beginning bias processing on {0} images".format(len(hdu_objects)))

    masks = []

    for hdu in hdu_objects:
        image_data = hdu.data
        center_quarter = extract_center_fraction_region(image_data, 0.25)
        center_quarter_median = numpy.median(center_quarter)

        corrected_image = image_data - center_quarter_median

        corrected_image_MAD = astropy.stats.median_absolute_deviation(corrected_image)

        sigma_range_start = corrected_image_MAD-(sigma_threshold * mad_constant * corrected_image_MAD)
        sigma_range_end = corrected_image_MAD+(sigma_threshold * mad_constant * corrected_image_MAD)
        # sigma_range_end = -sigma_range_start

        stddev_filtered_image = (corrected_image <= sigma_range_start) | \
                                (corrected_image >= sigma_range_end)

        masks.append(stddev_filtered_image)

    filtered_mask = apply_frequency_thresholding_on_masked_arrays(masks, 0.30)

    my_logger.debug("Completed bias processing")
    return filtered_mask


def apply_darks_processing(hdu_objects, dark_current_threshold=35):
    r"""**Algorithm**

    1. Locate the overscan region of each image, and divide the entire image by the median of the overscan region.

    2. Then, divide the exposure time of each image, to obtain the dark current, in units of electrons/second.

    3. Filter any images that have a dark current of more than 35 electrons/second.

    :param hdu_objects: A list of Header Data Unit objects.
    :param dark_current_threshold: The minimum dark current (in Electrons/second) you'd like to allow in your image.

    :return: A list of boolean arrays

    """

    my_logger.debug("Beginning darks processing on {0} images".format(
        len(hdu_objects)))

    masks = []

    for hdu in hdu_objects:
        bias_section_header_string = hdu.header['BIASSEC']
        exposure_time = float(hdu.header['EXPTIME'])
        image_data = hdu.data

        overscan_region_coordinates = get_slices_from_image_section(bias_section_header_string)

        overscan_region_median = numpy.median(image_data[overscan_region_coordinates])

        image_data = numpy.delete(image_data, overscan_region_coordinates[1], 1)
        image_data -= overscan_region_median
        image_data /= exposure_time
        masked_image = image_data > dark_current_threshold

        masks.append(masked_image)

    filtered_masks = apply_frequency_thresholding_on_masked_arrays(masks, 0.30)

    my_logger.debug("Completed darks processing.")
    return filtered_masks


def apply_flats_processing(hdu_objects, sigma_threshold=7):
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

    Where for a normal distribution, :math:`k\approx 1.4826`


    :param hdu_objects: A list of HDU objects

    """

    my_logger.debug("Beginning flats processing on {0} images".format(
        len(hdu_objects)))

    corrected_images_list = []

    for hdu in hdu_objects:
        image_data = hdu.data
        center_quarter = extract_center_fraction_region(hdu.data, 0.25)
        center_quarter_median = numpy.median(center_quarter)
        corrected_image = image_data / center_quarter_median
        corrected_images_list.append(corrected_image)

    # Create a 3D array out of all the corrected images, where (x,y) plane is the original image and the z-axis is what
    # connects the images
    corrected_images_cube = numpy.dstack(tuple(corrected_images_list))

    # Take the standard deviation of each pixel across all images
    # remember that axes are 0-indexed
    std_deviations_array = numpy.std(corrected_images_cube, axis=2)

    mad = astropy.stats.median_absolute_deviation(std_deviations_array)

    # sigma_range_start = -sigma_threshold * mad_constant * corrected_image_MAD

    # once you have the MAD, mask any values outside the range of the sigma threshold
    # range_start, range_end = mad - ((mad_constant * mad ) * sigma_threshold), \
    #                          mad + ((mad_constant * mad) * sigma_threshold)

    range_start, range_end = mad - ((mad_constant * mad ) * sigma_threshold), \
                             mad + ((mad_constant * mad) * sigma_threshold)

    filtered_array = numpy.ma.masked_outside(std_deviations_array, range_start,
                                             range_end)

    my_logger.debug("Completed flats processing")
    return numpy.ma.getmaskarray(filtered_array)


# ------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------
def extract_center_fraction_region(image, inner_edge_width):
    r"""Extract and return the center fraction of an image

    :param inner_edge_width: Size of inner edge as fraction of total image size
    :param image: image data
    :return: The center fraction of the original image
    :rtype: numpy.ndarray

    """
    ny, nx = image.shape

    inner_nx = round(nx * inner_edge_width)
    inner_ny = round(ny * inner_edge_width)
    return image[inner_ny: -inner_ny, inner_nx: -inner_nx]

def combine_image_masks(masks_by_type):
    """
    Once you have all N masks, where each mask corresponds to each type (i.e. you have 1 bias mask, 1 flat mask, etc.)
    then you can combine them to only extract

    :param masks_by_type: A list of image masks (boolean arrays)
    :return: numpy.ndarray
    """

    return numpy.logical_or.reduce(masks_by_type).astype(numpy.uint8)

def get_slices_from_image_section(image_section_string):
    """
    Borrowed from BANZAI. Convert FITS header image section value to tuple of slices.

    Example:  '[3100:3135,1:2048]' --> (slice(0, 2048, 1), slice(3099, 3135, 1))
    Note:
    FITS Header image sections are 1-based and indexed by [column, row]
    Numpy arrays are zero-based and indexed by [row, column]

    :param header_string: An image section string in the form "[x1:x2, y1:y2]"
    :return: Row-indexed tuple of slices, (row_slice, col_slice)
    """

    # Strip off the brackets and split the coordinates
    pixel_sections = image_section_string[1:-1].split(',')
    x_slice = split_slice(pixel_sections[0])
    y_slice = split_slice(pixel_sections[1])
    pixel_slices = (y_slice, x_slice)
    return pixel_slices

def split_slice(pixel_section):
    """
    Borrowed from BANZAI. Convert FITS header pixel section to Numpy-friendly
    slice.

    Example: '3100:3135' --> slice(3099, 3135, 1)
    """
    pixels = pixel_section.split(':')
    if int(pixels[1]) > int(pixels[0]):
        pixel_slice = slice(int(pixels[0]) - 1, int(pixels[1]), 1)
    else:
        if int(pixels[1]) == 1:
            pixel_slice = slice(int(pixels[0]) - 1, None, -1)
        else:
            pixel_slice = slice(int(pixels[0]) - 1, int(pixels[1]) - 2, -1)
    return pixel_slice

def apply_frequency_thresholding_on_masked_arrays(list_of_arrays,
                                                  frequency_threshold):
    """
    Takes a list of arrays and removes any values in the array don't appear more than a specified number of times.


    :param list_of_arrays: A list of numpy arrays with all the same shape
    :param frequency_threshold: A minimum frequency in the range (0, 1]
    :return: An array where only the values that appear more than the frequency threshold are preserved
    :rtype: numpy.ndarray
    """

    my_logger.debug("Applying frequency thresholding on {0} arrays.".format(
        len(list_of_arrays)))

    # if all arrays have the same shape as the first one then all arrays have the same shape

    if not (all(
        [array.shape == list_of_arrays[0].shape for array in list_of_arrays])):
        my_logger.warn(
            "The arrays you're attempting to apply frequency thresholding to do not all have the same size."
        )

    # sum the boolean arrays all together, the value at each (row, col) tells you the number of times a bad pixel
    # was detected there
    all_arrays_sum = sum(list_of_arrays)

    frequency_filtered_array = all_arrays_sum >= (frequency_threshold * len(list_of_arrays))

    return frequency_filtered_array
