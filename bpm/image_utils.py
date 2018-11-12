import numpy as np
import astropy.stats

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

    return np.uint8(outlier_mask)


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
