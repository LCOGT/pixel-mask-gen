import astropy.io.fits as fits
import os
import datetime
import glob
import logging
import lcogt_logging
import argparse
import bpm.image_processing as image_processing
import numpy as np

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

        logger.info("Found {num_frames} calibration frames to process".format(num_frames = len(calibration_frames)))

        dark_mask = image_processing.process_dark_frames(get_frames_of_type(calibration_frames, 'DARK'), args.dark_current_threshold)
        bias_mask = image_processing.process_bias_frames(get_frames_of_type(calibration_frames, 'BIAS'), args.bias_sigma_threshold)

        flats_sorted = sort_flats_by_filter((get_frames_of_type(calibration_frames, 'FLAT')))
        flat_masks = [image_processing.process_flat_frames(flats_sorted[filter], args.flat_sigma_threshold) for filter in flats_sorted.keys()]

        flat_masks.extend([bias_mask, dark_mask])
        combined_mask = np.sum(np.dstack(flat_masks), axis=2) > 0

        today_date = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S")
        instrument_code = calibration_frames[0].header['INSTRUME']
        header_info = {'OBSTYPE': 'BPM'}

        output_filename = os.path.join(args.output_directory, "bpm-{instrument}-{today}.fits".format(today=today_date,
                                                                                                     instrument=instrument_code))
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


def sort_flats_by_filter(frames):
    """
    Given a set of flat frames, sort them by filter and return a dict of the
    form:
    {filter_type:[frames_with_filter]}
    """
    filters = set([frame.header['FILTER'] for frame in frames])
    return {filter: [frame for frame in frames if frame.header['FILTER'] == filter]
                    for filter in filters}
