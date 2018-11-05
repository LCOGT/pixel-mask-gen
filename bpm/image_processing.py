import numpy as np
import astropy.stats

def process_bias_frames(bias_frames, mask_threshold=11):
    corrected_frames = []

    for frame in bias_frames:
        image_data = frame.data
        overscan_section = get_slices_from_header_section(frame.header['BIASSEC'])

        image_data -= np.median(image_data[overscan_section])

        corrected_frames.append(image_data)

    return mask_outliers(np.dstack(corrected_frames), mask_threshold)


def process_dark_frames(dark_frames, dark_current_threshold=35, mask_threshold=10):
    corrected_frames = []

    for frame in dark_frames:
        image_data = frame.data
        overscan_section = get_slices_from_header_section(frame.header['BIASSEC'])

        overscan_median = np.median(image_data[overscan_section])
        image_data -= overscan_median
        image_data /= float(frame.header['EXPTIME'])

        corrected_frames.append(image_data)

    dark_current_mask = np.mean(np.dstack(corrected_frames), axis=2) > dark_current_threshold
    outlier_mask = mask_outliers(np.dstack(corrected_frames), mask_threshold)

    return np.logical_or(dark_current_mask, outlier_mask)


def process_flat_frames(flat_frames, mask_threshold=11):
    filters = set([frame.header['FILTER'] for frame in flat_frames])

    if len(filters) != 1:
        raise ValueError("Flat frames are not of the same filter")

    corrected_frames = []

    for frame in flat_frames:
        image_data = frame.data
        overscan_section = get_slices_from_header_section(frame.header['BIASSEC'])
        trimsec_section = get_slices_from_header_section(frame.header['TRIMSEC'])

        overscan_median = np.median(image_data[overscan_section])
        image_median = np.median(image_data[trimsec_section])

        image_data -= overscan_median
        image_data /= image_median

        corrected_frames.append(image_data)

    return mask_outliers(np.dstack(corrected_frames), mask_threshold)

# ------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------

def mask_outliers(stacked_frames, mask_threshold=10):
    """
    Mask pixels outside a specified number of standard deviations

    Generate MAD for each pixel - 2D array - pixel_mads
    Generate standard deviation of all pixel MADs - Scalar - std_all_pixels
    Generate median of all pixels MADs - Scalar - median_all_pixels

    Flag any pixels whose std value is outside the median +/- mask_threshold * std_all_pixels

    :param stacked_frames: stack of corrected frames
    :param mask_threshold: standard deviation threshold
    """
    pixel_mads = astropy.stats.median_absolute_deviation(stacked_frames, axis=2)
    std_all_pixels = np.std(pixel_mads)
    median_all_pixels = np.median(pixel_mads)

    outlier_mask = np.logical_or(pixel_mads < median_all_pixels - (mask_threshold * std_all_pixels),
                                 pixel_mads > median_all_pixels + (mask_threshold * std_all_pixels))

    return outlier_mask


def get_slices_from_header_section(header_section_string):
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
    pixel_sections = header_section_string[1:-1].split(',')
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
