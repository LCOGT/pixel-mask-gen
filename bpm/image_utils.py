import numpy as np
import astropy.stats
import astropy.io.fits as fits

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


def get_extensions_by_name(fits_hdulist, name):
    """
    Get a list of the science extensions from a multi-extension fits file (HDU list)
    Parameters
    ----------
    fits_hdulist: HDUList
                  input fits HDUList to search for SCI extensions
    name: str
          Extension name to collect, e.g. SCI
    Returns
    -------
    HDUList: an HDUList object with only the SCI extensions
    """
    # The following of using False is just an awful convention and will probably be
    # deprecated at some point
    extension_info = fits_hdulist.info(False)
    return fits.HDUList([fits_hdulist[ext[0]] for ext in extension_info if ext[1] == name])


def apply_header_value_to_all_extensions(frames, header_keyword):
    """
    Apply a header value from an image's PrimaryHDU to its
    extensions.
    """
    for frame in frames:
        header_value = frame[0].header[header_keyword]
        for extension_num in range(1, len(frame)):
            frame[extension_num].header[header_keyword] = header_value


def stack_frames_by_extver(frames):
    """
    Given a 2D list of SCI extensions, flatten it and
    stack frames according to their EXTVER header keyword

    {header['EXTVER']: [frames_of_same_extver]}
    """
    frames_flattened = [fits.HDUList(extension) for extension in frames for extension in extension]
    return sort_frames_by_header_values(frames_flattened, 'EXTVER')



def sort_frames_by_header_values(frames, header_keyword):
    """
    Given a set of frames and a header keyword, sort the frames by the corresponding
    header values into a form:
    {header_value:[frames_with_header_value]}
    """
    header_values = set([frame[0].header[header_keyword] for frame in frames])
    return {value: [frame for frame in frames if frame[0].header[header_keyword] == value]
                    for value in header_values}
