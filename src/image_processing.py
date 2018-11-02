import numpy as np
import astropy.stats

def process_bias_frames(bias_frames, mad_threshold=8):
    corrected_frames = []

    for frame in bias_frames:
        image_data = frame.data
        overscan_coords = get_slices_from_image_section(frame.header['BIASSEC'])

        overscan_median = np.median(image_data[overscan_coords])
        image_data -= overscan_median

        corrected_frames.append(image_data)

    return mask_outliers(np.dstack(corrected_frames), mad_threshold)

def process_dark_frames(dark_frames, dark_current_threshold=35):
    corrected_frames = []

    for frame in dark_frames:
        image_data = frame.data
        overscan_region_coordinates = get_slices_from_image_section(frame.header['BIASSEC'])

        overscan_median = np.median(image_data[overscan_region_coordinates])
        image_data -= overscan_median
        image_data /= float(frame.header['EXPTIME'])

        corrected_frames.append(image_data)

    return np.mean(np.dstack(corrected_frames), axis=2) > dark_current_threshold

def process_flat_frames(flat_frames, mad_threshold=7):
    filters = set([frame.header['FILTER'] for frame in flat_frames])

    if (len(filters) != 1):
        raise ValueError("Flat frames are not of the same filter")

    corrected_frames = []

    for frame in flat_frames:
        image_data = frame.data
        overscan_coords = get_slices_from_image_section(frame.header['BIASSEC'])
        trimsec_coords = get_slices_from_image_section(frame.header['TRIMSEC'])

        overscan_median = np.median(image_data[overscan_coords])
        image_median = np.median(image_data[trimsec_coords])

        image_data -= overscan_median
        image_data /= image_median

        corrected_frames.append(image_data)

    return mask_outliers(np.dstack(corrected_frames), mad_threshold)

# ------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------

def mask_outliers(stacked_frames, num_mads=7):
    """
    Mask pixels outside a specified number of median absolute deviations

    Generate MAD for each pixel - 2D array - mad_array
    Generate MAD of all pixel MADs - Scalar - mad_of_data
    Generate median of all pixels MADs - Scalar - median_of_data

    Flag any pixels whose MAD is outside the median +/- num_mads * mad_of_data
    """
    mad_array = astropy.stats.median_absolute_deviation(stacked_frames, axis=2)
    mad_of_data = astropy.stats.median_absolute_deviation(mad_array)
    median_of_data = np.median(mad_array)

    outlier_mask = np.logical_or(mad_array < median_of_data - (num_mads * mad_of_data),
                                 mad_array > median_of_data + (num_mads * mad_of_data))
    return outlier_mask


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
