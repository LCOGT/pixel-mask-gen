import astropy.io.fits as fits
import os
import datetime
import glob
import logging
import lcogt_logging
import argparse
import bpm.image_processing as image_processing
import numpy as np
import bpm.image_utils as image_utils

logger = logging.getLogger('lco-bpm-maker')

def setup_logging(log_level):
    logger.setLevel(log_level)
    handler = logging.StreamHandler()
    handler.setLevel(log_level)
    handler.setFormatter(lcogt_logging.LCOGTFormatter())
    logger.addHandler(handler)


def parse_args():
    parser = argparse.ArgumentParser(description='Create a bad pixel mask from a set of calibration frames.')
    parser.add_argument('input_directory', help='Input directory of calibration images')
    parser.add_argument('output_directory', help='Output directory for bad pixel mask')
    parser.add_argument('--log-level', dest='log_level', default='INFO', help='Logging level to be displayed',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])
    parser.add_argument('--dark-current-threshold', dest='dark_current_threshold',
                        help='Threshold for pixel dark current when flagging bad pixels in dark frames. Pixels above this will be flagged. Default = 35 [electrons/second]',
                        default=35)
    parser.add_argument('--flat-sigma-threshold', dest='flat_sigma_threshold',
                        help='Number of standard deviations from the median of the combined flat image for a pixel to be flagged as bad. Default = 10',
                        default=10)
    parser.add_argument('--bias-sigma-threshold', dest='bias_sigma_threshold',
                        help='Number of standard deviations from the median of the combined bias image for a pixel to be flagged as bad. Default = 10',
                        default=10)

    args = parser.parse_args()

    return args


def generate_bpm():
    args = parse_args()
    setup_logging(getattr(logging, args.log_level))

    calibration_frames = get_calibration_frames(os.path.normpath(args.input_directory) + '/*.fits')

    if calibration_frames:
        if len(set([frame[0].header['INSTRUME'] for frame in frames])) != 1:
            raise RuntimeError("Got calibration frames from more than one camera. Aborting.")

        #separate out single and multi-extension FITS files
        single_extension_frames = [frame[0] for frame in calibration_frames if len(frame) == 1]
        multi_extension_frames = [frame for frame in calibration_frames if frame not in single_extension_frames]

        #Create BPMs!
        process_multi_extension_frames(multi_extension_frames,
                                       int(args.dark_current_threshold),
                                       int(args.flat_sigma_threshold),
                                       ing(args.bias_sigma_threshold))

        process_single_extension_frames(single_extension_frames,
                                        int(args.dark_current_threshold),
                                        int(args.flat_sigma_threshold),
                                        ing(args.bias_sigma_threshold))
    else:
        raise RuntimeError("No calibration frames could be found. Check that the directory contains calibration frames.")


def process_multi_extension_frames(frames,
                                   dark_current_threshold,
                                   flat_sigma_threshold,
                                   bias_sigma_threshold):
    if len(frames) == 0:
        return

    camera_has_no_overscan = True
    try:
        bias_section = image_utils.get_slices_from_header_section(frames[0][1].header['BIASSEC'])
        camera_has_no_overscan = False
    except:
        logger.warn("Couldn't parse BIASSEC keyword. Using bias frames to determine camera bias level.")

    #Sort frames by binning
    frames_sorted_by_binning = image_utils.sort_frames_by_header_values(frames, 'CCDSUM')

    logger.info("Beginning processing on {num_frames} calibration frames".format(num_frames = len(frames)))
    for binning in frames_sorted_by_binning.keys():
        frames_sorted = frames_sorted_by_binning[binning]
        #TODO: make sure this actually works
        bias_level  = image_processing.get_bias_level_from_frames(get_frames_of_type(frames_sorted_by_binning[binning], 'BIAS')) if camera_has_no_overscan else None

        #Update SCI extensions with required header values for image processing methods
        image_utils.apply_header_value_to_all_extensions(get_frames_of_type(frames_sorted, 'DARK'), 'EXPTIME')
        image_utils.apply_header_value_to_all_extensions(get_frames_of_type(frames_sorted, 'FLAT'), 'FILTER')

        #Retrieve all SCI extensions for calibration frames
        dark_sci_extensions = [get_extensions_by_name(frame, 'SCI') for frame in get_frames_of_type(frames_sorted, 'DARK')]
        bias_sci_extensions = [get_extensions_by_name(frame, 'SCI') for frame in get_frames_of_type(frames_sorted, 'BIAS')]
        flat_sci_extensions = [get_extensions_by_name(frame, 'SCI') for frame in get_frames_of_type(frames_sorted, 'FLAT')]

        #For flats, separate SCI extensions by filter
        flat_sci_extensions_sorted = image_utils.sort_frames_by_header_values(flat_sci_extensions, 'FILTER')

        #for every filter type, further sort flat SCI extensions by EXTVER keyword
        flat_stacks = []
        for filter in flat_sci_extensions_sorted.keys():
            flat_stacks.append(stack_frames_by_extver(flat_sci_extensions_sorted[filter]))

        #Stack bias and dark frames according to EXTVER
        dark_frames_stacked = image_utils.stack_frames_by_extver(dark_sci_extensions)
        bias_frames_stacked = image_utils.stack_frames_by_extver(bias_sci_extensions)

        #Now that we have all images stacked, create bpms!
        dark_bpms = {extver_key: image_processing.process_dark_frames([frame[0] for frame in dark_frames_stacked[extver_key]],
                                                                      dark_current_threshold,
                                                                      bias_level)
                     for extver_key in dark_frames_stacked.keys()}

        bias_bpms = {extver_key: image_processing.process_bias_frames([frame[0] for frame in bias_frames_stacked[extver_key]],
                                                                      bias_sigma_threshold)
                     for extver_key in dark_frames_stacked.keys()}}

        flat_bpms = []
            for filter in flat_stacks:
                flat_bpms_for_filter = {extver_key: image_processing.process_flat_frames([frame[0] for frame in filter[extver_key]],
                                                                                         flat_sigma_threshold)
                                        for extver_key in filter.keys()}
                flat_bpms.append(flat_bpms_for_filter)


def process_single_extension_frames(frames,
                                    dark_current_threshold,
                                    flat_sigma_threshold,
                                    bias_sigma_threshold):
    if len(frames) == 0:
        return

    #check if camera has an overscan region
    camera_has_no_overscan = True
    try:
        bias_section = image_utils.get_slices_from_header_section(frames[0][0].header['BIASSEC'])
        camera_has_no_overscan = False
    except:
        logger.warn("Couldn't parse BIASSEC keyword. Using bias frames to determine camera bias level.")

    frames_sorted_by_binning = image_utils.sort_frames_by_header_values(frames, 'CCDSUM')

    logger.info("Beginning processing on {num_frames} calibration frames".format(num_frames = len(frames)))
    for binning in frames_sorted_by_binning.keys():
        frames_sorted = frames_sorted_by_binning[binning]

        bias_level  = image_processing.get_bias_level_from_frames(get_frames_of_type(frames_sorted, 'BIAS')) if camera_has_no_overscan else None

        dark_mask = image_processing.process_dark_frames(get_frames_of_type(frames_sorted, 'DARK'), int(dark_current_threshold), bias_level)
        bias_mask = image_processing.process_bias_frames(get_frames_of_type(frames_sorted, 'BIAS'), int(bias_sigma_threshold))

        flats_sorted = image_utils.sort_frames_by_header_values((get_frames_of_type(frames_sorted, 'FLAT')), 'FILTER')
        flat_masks = [image_processing.process_flat_frames(flats_sorted[filter], int(flat_sigma_threshold), bias_level) for filter in flats_sorted.keys()]

        flat_masks.extend([bias_mask, dark_mask])
        combined_mask = np.sum(np.dstack(flat_masks), axis=2) > 0

        today_date = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H%M%S")
        instrument_code = frames_sorted_by_binning[binning][0].header['INSTRUME']

        header_info = {'OBSTYPE': 'BPM',
                       'DAY-OBS': frames_sorted_by_binning[binning][0].header['DAY-OBS'],
                       'CCDSUM': binning,
                       'SITEID': frames_sorted_by_binning[binning][0].header['SITEID'],
                       'INSTRUME': instrument_code,
                       'DATE-OBS': datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]}

        output_filename = os.path.join(args.output_directory, "bpm-{instrument}-{bin_type}-{today}.fits".format(instrument=instrument_code,
                                                                                                                today=today_date,
                                                                                                                bin_type=binning.replace(" ", "x")))
        fits.writeto(filename=output_filename,
                     data=combined_mask.astype(np.uint8),
                     header=fits.Header(header_info),
                     checksum=True)

        logger.info("Finished processing. BPM written to {file_path}".format(file_path=output_filename))


def get_calibration_frames(path_to_frames, calibration_types=['d00', 'f00', 'b00']):
    """
    Given a directory of fits files, return a list of calibration frames
    """
    frames = glob.glob(path_to_frames)
    frames = [fits.open(frame, mode='readonly') for frame in frames if any(obs_type in frame for obs_type in calibration_types)]
    return frames


def get_frames_of_type(frames, observation_type):
    """
    Takes in a list of frames, and returns frames which match the observation
    type provided.
    """
    return [frame for frame in frames if observation_type in frame[0].header['OBSTYPE']]
