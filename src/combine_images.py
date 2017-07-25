
from sys import argv
#from os import listdir
from glob import glob
from astropy.io import fits
import numpy as np
import fitsutils

#####
# FUNCTIONS

def create_zeroed_image(image):

    hdu = fits.open("{0}".format(image))
    Y_axis = hdu[0].header['NAXIS2']
    X_axis = hdu[0].header['NAXIS1']
    hdu.close()

    array_of_zeros = np.zeros((Y_axis,X_axis))

    return array_of_zeros


def reassign_obstype(image):

    hdu = fits.open("{0}".format(image))
    hdu[0].header['OBSTYPE'] = 'BPM'

    imageheader = hdu[0].header
    imagedata = hdu[0].data

    hdu.close

    return imageheader, imagedata

    
#####
# COMMAND LINE
if __name__ == '__main__':

    camera = argv[1]

    bpmlist = glob("{0}/{0}_bpm.*.fits".format(camera))

    for bpm in bpmlist:

        if bpm == bpmlist[0]:
        
            final_bpm = create_zeroed_image(bpm) # initialize array for BPM

        retrieved_header, retrieved_bpm = reassign_obstype(bpm)

        final_bpm = final_bpm + retrieved_bpm

    #final_image[513,758] = 1 # extra bad pix for kb88
    final_bpm[final_bpm > 0] = 1
    final_bpm[final_bpm < 0] = 0

    iexec = fitsutils.outputFITS(final_bpm,retrieved_header,'combined_bpm.fits')
