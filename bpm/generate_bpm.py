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
        if len(set([frame[0].header['INSTRUME'] for frame in calibration_frames])) != 1:
            raise RuntimeError("Got calibration frames from more than one camera. Aborting.")

        #separate out single and multi-extension FITS files
        multi_extension_frames = [frame for frame in calibration_frames if len(frame) > 1]
        single_extension_frames = [frame[0] for frame in calibration_frames if frame not in multi_extension_frames]

        #Create BPMs!
        process_multi_extension_frames(multi_extension_frames,
                                       int(args.dark_current_threshold),
                                       int(args.flat_sigma_threshold),
                                       int(args.bias_sigma_threshold),
                                       args.output_directory)

        process_single_extension_frames(single_extension_frames,
                                        int(args.dark_current_threshold),
                                        int(args.flat_sigma_threshold),
                                        int(args.bias_sigma_threshold),
                                        args.output_directory)
    else:
        raise RuntimeError("No calibration frames could be found. Check that the directory contains calibration frames.")


def process_multi_extension_frames(frames,
                                   dark_current_threshold,
                                   flat_sigma_threshold,
                                   bias_sigma_threshold,
                                   output_directory):
    """
    Process a set of multi-extension FITS frames into bad pixel mask(s).
    """
    if len(frames) == 0:
        return

    #determine number of amplifiers from image
    n_amplifiers = len(image_utils.get_extensions_by_name(frames[0], 'SCI'))

    camera_has_no_overscan = True
    try:
        bias_section = image_utils.get_slices_from_header_section(frames[0][1].header['BIASSEC'])
        camera_has_no_overscan = False
    except:
        logger.warn("Couldn't parse BIASSEC keyword. Using bias frames to determine camera bias level.")

    #Sort frames by binning
    frames_sorted_by_binning = image_utils.sort_frames_by_header_values(frames, 'CCDSUM')

    logger.info("Beginning processing on {num_frames} calibration frames".format(num_frames = len(frames)))

    #For each binning scheme, create a separate BPM
    for binning in frames_sorted_by_binning.keys():
        frames_sorted = frames_sorted_by_binning[binning]

        #Update SCI extensions with required header values for image processing methods
        header_keywords_to_update = ['OBSTYPE', 'EXPTIME', 'FILTER']

        for keyword in header_keywords_to_update:
            image_utils.apply_header_value_to_all_extensions(multi_extension_frames, keyword)

        #Treat each amplifier as a separate image
        masks = []
        for amplifier in range(0, n_amplifiers):
            amplifier_frames = image_utils.get_sci_extensions_from_amplifier(frames_sorted, 0)

            combined_mask = create_final_mask(frames_sorted,
                                              dark_current_threshold,
                                              flat_sigma_threshold,
                                              bias_sigma_threshold,
                                              camera_has_no_overscan)

            masks.append(combined_mask)

        header = fits.Header({'OBSTYPE': 'BPM',
                              'DAY-OBS': frames_sorted_by_binning[binning][0].header['DAY-OBS'],
                              'CCDSUM': binning,
                              'SITEID': frames_sorted_by_binning[binning][0].header['SITEID'],
                              'INSTRUME': frames_sorted_by_binning[binning][0].header['INSTRUME'],
                              'DATE-OBS': datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]})

        write_bpm_to_file(masks, output_directory, header)


def process_single_extension_frames(frames,
                                    dark_current_threshold,
                                    flat_sigma_threshold,
                                    bias_sigma_threshold,
                                    output_directory):
    """
    Process a set of single-extension FITS frames into bad pixel mask(s).
    """
    if len(frames) == 0:
        return

    #check if camera has an overscan region
    camera_has_no_overscan = True
    try:
        bias_section = image_utils.get_slices_from_header_section(frames[0].header['BIASSEC'])
        camera_has_no_overscan = False
    except:
        logger.warn("Couldn't parse BIASSEC keyword. Using bias frames to determine camera bias level.")

    frames_sorted_by_binning = image_utils.sort_frames_by_header_values(frames, 'CCDSUM')

    logger.info("Beginning processing on {num_frames} calibration frames".format(num_frames = len(frames)))
    for binning in frames_sorted_by_binning.keys():
        frames_sorted = frames_sorted_by_binning[binning]

        combined_mask = create_final_mask(frames_sorted,
                                          dark_current_threshold,
                                          flat_sigma_threshold,
                                          bias_sigma_threshold,
                                          camera_has_no_overscan)

        header = fits.Header({'OBSTYPE': 'BPM',
                              'DAY-OBS': frames_sorted_by_binning[binning][0].header['DAY-OBS'],
                              'CCDSUM': binning,
                              'SITEID': frames_sorted_by_binning[binning][0].header['SITEID'],
                              'INSTRUME': frames_sorted_by_binning[binning][0].header['INSTRUME'],
                              'DATE-OBS': datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]})

        write_bpm_to_file(combined_mask, output_directory, header)


def write_bpm_to_file(combined_masks, output_directory, header):
    """
    Write output BPM to file.
    """

    today_date = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H%M%S")
    output_filename = os.path.join(output_directory, "bpm-{instrument}-{bin_type}-{today}.fits".format(instrument=header['INSTRUME'],
                                                                                                       today=today_date,
                                                                                                       bin_type=header['CCDSUM'].replace(" ", "x")))

    if combined_masks.ndim != 1:
        hdu_list = fits.HDUList([fits.ImageHDU(data=mask, header=header) for mask in combined_masks])
    else:
        fits.writeto(filename=output_filename,
                     data=combined_masks.astype(np.uint8),
                     header=header,
                     checksum=True)

    logger.info("Finished processing. BPM written to {file_path}".format(file_path=output_filename))


def create_final_mask(frames, dark_current_threshold, bias_sigma_threshold, flat_sigma_threshold, camera_has_no_overscan=True):
    """
    From a set of calibration frames, create a final bad pixel mask.
    """
    bias_level  = image_processing.get_bias_level_from_frames(get_frames_of_type(frames, 'BIAS')) if camera_has_no_overscan else None

    dark_mask = image_processing.process_dark_frames(get_frames_of_type(frames, 'DARK'), int(dark_current_threshold), bias_level)
    bias_mask = image_processing.process_bias_frames(get_frames_of_type(frames, 'BIAS'), int(bias_sigma_threshold))

    flats_sorted = image_utils.sort_frames_by_header_values((get_frames_of_type(frames, 'FLAT')), 'FILTER')
    flat_masks = [image_processing.process_flat_frames(flats_sorted[filter], int(flat_sigma_threshold), bias_level) for filter in flats_sorted.keys()]

    flat_masks.extend([bias_mask, dark_mask])
    combined_mask = np.sum(np.dstack(flat_masks), axis=2) > 0

    return combined_mask


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
    return [frame for frame in frames if observation_type in frame.header['OBSTYPE']]
