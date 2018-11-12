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
        overscan_section = image_utils.get_slices_from_header_section(frame.header['BIASSEC'])
        trimsec_section = image_utils.get_slices_from_header_section(frame.header['TRIMSEC'])

        image_data -= np.median(image_data[overscan_section])

        corrected_frames.append(image_data[trimsec_section])

    return image_utils.mask_outliers(np.dstack(corrected_frames), mask_threshold)


def process_dark_frames(dark_frames, dark_current_threshold=35):

    logger.info("Processing {num_frames} dark frames".format(num_frames=len(dark_frames)))
    corrected_frames = []

    for frame in dark_frames:
        image_data = np.float32(frame.data)
        overscan_section = image_utils.get_slices_from_header_section(frame.header['BIASSEC'])
        trimsec_section = image_utils.get_slices_from_header_section(frame.header['TRIMSEC'])

        image_data -= np.median(image_data[overscan_section])
        image_data /= np.float32(frame.header['EXPTIME'])

        corrected_frames.append(image_data[trimsec_section])

    return np.uint8(np.median(np.dstack(corrected_frames), axis=2) > dark_current_threshold)


def process_flat_frames(flat_frames, mask_threshold=10):

    logger.info("Processing {num_frames} flat frames taken with filter: {filter}".format(num_frames=len(flat_frames),
                                                                                         filter=flat_frames[0].header['FILTER']))
    corrected_frames = []

    for frame in flat_frames:
        image_data = np.float32(frame.data)
        overscan_section = image_utils.get_slices_from_header_section(frame.header['BIASSEC'])
        trimsec_section = image_utils.get_slices_from_header_section(frame.header['TRIMSEC'])

        image_data -= np.median(image_data[overscan_section])
        image_data /= np.median(image_data[trimsec_section])

        corrected_frames.append(image_data[trimsec_section])

    return image_utils.mask_outliers(np.dstack(corrected_frames), mask_threshold)
