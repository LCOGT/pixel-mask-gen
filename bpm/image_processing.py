import numpy as np
import astropy.stats
import logging
import bpm.image_utils as image_utils

logger = logging.getLogger('lco-bpm-maker')

def process_bias_frames(bias_frames, mask_threshold=10):

    logger.info("Processing {num_frames} bias frames".format(num_frames=len(bias_frames)))
    corrected_frames = []

    for frame in bias_frames:
        image_data = np.float32(frame.data)
        trimsec_section = image_utils.get_slices_from_header_section(frame.header['TRIMSEC'])

        image_data -= np.median(image_data)

        corrected_frames.append(image_data[trimsec_section])

    return image_utils.mask_outliers(np.dstack(corrected_frames), mask_threshold)


def process_dark_frames(dark_frames, dark_current_threshold=35, bias_level=None):

    logger.info("Processing {num_frames} dark frames".format(num_frames=len(dark_frames)))
    corrected_frames = []

    for frame in dark_frames:
        image_data = np.float32(frame.data)
        if bias_level is None:
            overscan_section = image_utils.get_slices_from_header_section(frame.header['BIASSEC'])
            bias_level = np.median(image_data[overscan_section])

        trimsec_section = image_utils.get_slices_from_header_section(frame.header['TRIMSEC'])

        image_data -= bias_level
        image_data /= np.float32(frame.header['EXPTIME'])

        corrected_frames.append(image_data[trimsec_section])

    return np.uint8(np.median(np.dstack(corrected_frames), axis=2) > dark_current_threshold)


def process_flat_frames(flat_frames, mask_threshold=10, bias_level=None):

    logger.info("Processing {num_frames} flat frames taken with filter: {filter}".format(num_frames=len(flat_frames),
                                                                                         filter=flat_frames[0].header['FILTER']))
    corrected_frames = []

    for frame in flat_frames:
        image_data = np.float32(frame.data)
        if bias_level is None:
            overscan_section = image_utils.get_slices_from_header_section(frame.header['BIASSEC'])
            bias_level = np.median(image_data[overscan_section])

        trimsec_section = image_utils.get_slices_from_header_section(frame.header['TRIMSEC'])

        image_data -= bias_level
        image_data /= np.median(image_data[trimsec_section])

        corrected_frames.append(image_data[trimsec_section])

    return image_utils.mask_outliers(np.dstack(corrected_frames), mask_threshold)

def get_bias_level_from_frames(bias_frames):
    trimsec_section = image_utils.get_slices_from_header_section(bias_frames[0].header['TRIMSEC'])
    bias_level = np.median([np.median(frame.data[trimsec_section]) for frame in bias_frames])

    return bias_level
