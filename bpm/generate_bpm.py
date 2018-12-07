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
        if len(set([frame.header['INSTRUME'] for frame in calibration_frames])) != 1:
            raise RuntimeError("Got calibration frames from more than one camera. Aborting.")

        #check if camera has an overscan region
        camera_has_no_overscan = True
        try:
            bias_section = image_utils.get_slices_from_header_section(calibration_frames[0].header['BIASSEC'])
            camera_has_no_overscan = False
        except:
            logger.warn("Couldn't parse BIASSEC keyword. Using bias frames to determine camera bias level.")

        frames_sorted_by_binning = sort_frames_by_header_values(calibration_frames, 'CCDSUM')

        logger.info("Beginning processing on {num_frames} calibration frames".format(num_frames = len(calibration_frames)))
        for binning in frames_sorted_by_binning.keys():
            bias_level  = image_processing.get_bias_level_from_frames(get_frames_of_type(frames_sorted_by_binning[binning], 'BIAS')) if camera_has_no_overscan else None

            dark_mask = image_processing.process_dark_frames(get_frames_of_type(frames_sorted_by_binning[binning], 'DARK'), int(args.dark_current_threshold), bias_level)
            bias_mask = image_processing.process_bias_frames(get_frames_of_type(frames_sorted_by_binning[binning], 'BIAS'), int(args.bias_sigma_threshold))

            flats_sorted = sort_frames_by_header_values((get_frames_of_type(frames_sorted_by_binning[binning], 'FLAT')), 'FILTER')
            flat_masks = [image_processing.process_flat_frames(flats_sorted[filter], int(args.flat_sigma_threshold), bias_level) for filter in flats_sorted.keys()]

            flat_masks.extend([bias_mask, dark_mask])
            combined_mask = np.sum(np.dstack(flat_masks), axis=2) > 0

            today_date = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S")
            instrument_code = frames_sorted_by_binning[binning][0].header['INSTRUME']

            header_info = {'OBSTYPE': 'BPM',
                           'DAY-OBS': frames_sorted_by_binning[binning][0].header['DAY-OBS'],
                           'CCDSUM': binning,
                           'SITEID': frames_sorted_by_binning[binning][0].header['SITEID'],
                           'INSTRUME': instrument_code}

            output_filename = os.path.join(args.output_directory, "bpm-{instrument}-{bin_type}-{today}.fits".format(instrument=instrument_code,
                                                                                                                    today=today_date,
                                                                                                                    bin_type=binning.replace(" ", "x")))
            fits.writeto(filename=output_filename,
                         data=combined_mask.astype(np.uint8),
                         header=fits.Header(header_info),
                         checksum=True)

            logger.info("Finished processing. BPM written to {file_path}".format(file_path=output_filename))

    else:
        raise RuntimeError("No calibration frames could be found. Check that the directory contains calibration frames.")


def get_calibration_frames(path_to_frames, calibration_types=['d00', 'f00', 'b00']):
    """
    Given a directory of fits files, return a list of calibration frames
    """
    frames = glob.glob(path_to_frames)
    frames = [fits.open(frame, mode='readonly')[0] for frame in frames if any(obs_type in frame for obs_type in calibration_types)]
    return frames


def get_frames_of_type(frames, observation_type):
    """
    Takes in a list of frames, and returns frames which match the observation
    type provided.
    """
    return [frame for frame in frames if observation_type in frame.header['OBSTYPE']]

def sort_frames_by_header_values(frames, header_keyword):
    """
    Given a set of flat frames and a header keyword, sort the frames by the corresponding
    header values into a form:
    {header_value:[frames_with_header_value]}
    """
    header_values = set([frame.header[header_keyword] for frame in frames])
    return {value: [frame for frame in frames if frame.header[header_keyword] == value]
                    for value in header_values}
