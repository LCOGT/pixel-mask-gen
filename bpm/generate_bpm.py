import astropy.io.fits
import os
import sys
import datetime
import glob
import time
import logging
import lcogt_logging
import argparse
import bpm.image_processing as image_processing

logger = logging.getLogger('generate-bpm')

def setup_logging(log_level):
    logger.setLevel(log_level)
    handler = logging.StreamHandler()
    handler.setLevel(log_level)
    handler.setFormatter(lcogt_logging.LCOGTFormatter())
    logger.addHandler(handler)


def generate_bpm():

    parser = argparse.ArgumentParser(description='Create a bad pixel mask from a set of calibration frames.')
    parser.add_argument('--input-directory', dest='in', required=True, help='Input directory of calibration images')
    parser.add_argument('--output-directory', dest='out', required=True, help='Output directory for bad pixel mask')
    parser.add_argument('--log-level', dest='log_level', default='INFO', help='Logging level',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])

    args = parser.parse_args()
    setup_logging(getattr(logging, args.log_level))


    #
    # current_time = time.time()
    #
    # my_logger.info("Started main function")
    # absolute_source_file_path = os.path.abspath(source_file_path)
    #
    # hdu_for_images = get_image_hdus(
    #     absolute_source_file_path, suffix_list=["b00", "d00", "f00"])
    #
    # # This is a code smell, please fix this
    # bias_hdu_objects = sort_images_by_type('BIAS', hdu_for_images)
    # dark_hdu_objects = sort_images_by_type('DARK', hdu_for_images)
    # flat_hdu_objects = sort_images_by_type('SKYFLAT', hdu_for_images)
    #
    # bias_image_masks = image_processing.apply_bias_processing(bias_hdu_objects)
    #
    # dark_image_masks = image_processing.apply_darks_processing(
    #     dark_hdu_objects)
    #
    # flat_image_masks = image_processing.apply_flats_processing(
    #     flat_hdu_objects)
    #
    # # once you have the image masks, OR them together to get one mask
    # combined_image_mask = image_processing.combine_image_masks(
    #     [bias_image_masks, dark_image_masks, flat_image_masks])
    #
    # today = datetime.datetime.now()
    # format_strings = ["%Y-%m-%d", "-%H%M%S"]
    # today_date, today_time = tuple(
    #     [today.strftime(str_format) for str_format in format_strings])
    #
    # header_dict = {
    #     'OBSTYPE': 'BPM',
    #     'DATE': today.strftime("%Y-%m-%d"),
    #     'TIME': today.strftime("-%H%M%S")
    # }
    #
    # # cant directly write headers as dict object to astropy, must convert to header object first
    # output_filename = os.path.join(
    #     'output', "bpm-{0}-{1}.fits".format(today_date, today_time))
    # astropy.io.fits.writeto(
    #     filename=output_filename,
    #     data=combined_image_mask,
    #     header=astropy.io.fits.Header(header_dict))
    #
    # finished_time = time.time()
    #
    # time_diff = round(finished_time - current_time, 3)
    # my_logger.info("Time elapsed: {0}".format(time_diff))

def get_calibration_frames(path_to_frames, calibration_types=['BIAS', 'DARK', 'SKYFLAT', 'LAMPFLAT', 'DOMEFLAT']):
    """
    Given a directory of fits files, return a list of calibration frames
    """
    frames = glob.glob(path_to_frames)
    frames = [fits.open(frame, mode='readonly') for frame in frames if fits.getval(frame, 'OBSTYPE') in calibration_types]
    return frames

def get_frames_of_type(frames, observation_type):
    """
    Takes in a list of frames, and returns frames which match the observation
    type provided.
    """
    return [frame for frame in frames if frame.header['OBSTYPE'] == observation_type]

def sort_flats_by_filter(frames):
    """
    Given a set of flat frames, sort them by filter and return a dict of the
    form:
    {filter_type:[frames_with_filter]}
    """
    filters = set([frame['FILTER'] for frame in frames])
    return {filter: [frame for frame in frames if frame['FILTER'] == filter]
                    for filter in filters}
