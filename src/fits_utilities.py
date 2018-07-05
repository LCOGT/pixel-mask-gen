# External imports
import astropy.io.fits
import datetime
import os
import numpy
import logging

import pdb


def output_to_FITS(image_data, header_dict, filename, debug=False):
    """Generates a FITS v4 file from image data

    :param image_data: A numpy array, will be used for the primary header.
    :param header_dict: A dictionary used for the header data unit key/value pairs.
    :param filename: The destination file name for the FITS file.
    :return: None

    """
    # See: https://fits.gsfc.nasa.gov/fits_primer.html for info about FITS

    if image_data.size < 1:
        raise ValueError("No data to write to Primary HDU -- cannot read empty array.")

    else:
        new_hdu = astropy.io.fits.PrimaryHDU(image_data.astype(numpy.uint8))
        logging.info("Writing array of size {0} to file; array contains {1} bad pixels in mask.".format(image_data.shape,
                                                                                                       image_data.sum()))
        new_hdu_list = astropy.io.fits.HDUList([new_hdu])

    # Do you really need to set every header, or is it fine to just set the OBSTYPE header?
    new_hdu_list[0].header.set('OBSTYPE', 'BPM')

    todays_date = datetime.datetime.today().strftime("%Y%m%d")


    if debug == True:
        logging.info("Writing debugging information to /debug/")
        # For debugging purposes, write the mask and the header file into a text file
        #final_bpm_txtfile_path = os.path.join('debug', todays_date + "-" + filename + ".txt")
        final_bpm_txtfile_path = os.path.join("debug", todays_date + "_final-bmp.txt")

        with open(str(final_bpm_txtfile_path), 'w') as final_bpm_txtfile:
            final_bpm_txtfile.write(''.join(map(str, [coords for coords in image_data])))

            if not header_dict:
                final_bpm_txtfile.write('==========')
                final_bpm_txtfile.write('HEADERS\n')
                header_string = ''
                for key, value in header_dict.items():
                    final_bpm_txtfile.write("key: {0}, value: {1}\n".format(key, value))
                final_bpm_txtfile.write(header_string)

            final_bpm_txtfile.close()

    final_filename = filename + "-" + todays_date + "_bpm.fits"

    logging.info("Preparing to write data to file: {0}".format(final_filename))
    new_hdu_list.writeto(final_filename , overwrite=True, output_verify='exception',checksum=True)

    new_hdu_list.close()

    # after closing the file, ensure that the file it wrote to was not empty.
    if os.stat(filename).st_size < 10:
        raise ValueError('It appears the FITS file given from the data could not be written to.')


def read_individual_fits_file(filename):
    """

    :param filename: The absolute filename of a FITS image


    :returns image_data: A numpy array representing the data in the image
    :returns image_header_info: A dictionary object representing the header information for the last image. The \
                important header information will not change from image to image usually
    :returns image_shape: A tuple in the form (x,y) where x denotes the number of rows and y denotes the number of columns

    """

    if len(filename) < 1:
        raise ValueError('Filename is blank.')

    logging.info("Reading FITS file: {0}".format(filename))

    image_file = astropy.io.fits.open(filename)
    image_data = image_file[0].data
    image_header_info = image_file[0].header
    image_shape = image_data.shape

    image_file.close()

    return image_data, image_header_info, image_shape
