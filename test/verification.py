"""
This script is used to parse the bad pixel mask from BANZAI-processed images that were\
downloaded from the LCO archive at: https://lco.global/observatory/data/archive

The purpose of this script is to be able to test your bad pixel mask against a known-good bad\
pixel mask during the development process.

After Banzai processes images, it stores the bad pixel mask for them in the 2nd extension of the\
FITS header data unit.
"""

import astropy.io.fits
import pdb
import os
import numpy


def read_bpms_from_files(directory):
    bpms = []
    print('Searching for FITS files in directory.')
    for index, filename in enumerate(os.listdir(".")):
        if filename.endswith('.fits'):
            print("index:{0}, filename:{1}".format(index, filename))
            image_file = astropy.io.fits.open(filename)
            # BANZAI stores the BPM in the 2nd extension of the image
            # https://lco.global/documentation/archive-documentation/
            bpm_data = image_file[1].data

            bpms.append(bpm_data)

    if len(bpms) < 1:
        raise ValueError("No FITS images found in directory: {0}".format(os.getcwd()))
    return bpms

def verify_bpms_are_equal(bpms):
    print('Checking that all bad pixels masks in the 2nd extension of the FITS images are equal.')

    if len(bpms) == 2:
        return numpy.array_equal(bpms[0], bpms[1])

    else:
        # Check if all possible pairs of numpy arrays are equal
        all([numpy.array_equal(bpms[i], bpms[j]) for i in range(len(bpms)) for j in range(i+1, len(bpms))])

if __name__ == '__main__':
    # The directory where the BANZAI-reduced images are stored
    DEFAULT_DIRECTORY = 'sample_images/processed/'
    bpms_from_directory = read_bpms_from_files()
    if verify_bpms_are_equal(bpms_from_directory):
        '''
        print("Writing bpm to file.")
        with open('original_bpm.txt', 'w') as file_handler:
            file_handler.write(''.join(map(str, [coords for coords in bpms_from_directory[0]])))
    file_handler.close()
        '''
        print('Converting BPM to FITS file.')

        new_hdu = astropy.io.fits.PrimaryHDU(bpms_from_directory[0].astype(numpy.uint8))
        new_hdu_list = astropy.io.fits.HDUList([new_hdu])

        new_hdu_list[0].header.set('OBSTYPE', 'BPM')

        filename = 'sample_bpm.fits'
        new_hdu_list.writeto(filename, overwrite=True, output_verify='exception', checksum=True)
        new_hdu_list.close()

    else:
        print('Not all FITS files had the same bad pixel in the 2nd extension mask.')



