'''
    CREATE_BPM.PY

    A script to create a bad pixel mask (BPM) from a set of calibration frames.
'''

#####
# IMPORTED FUNCTIONS
from astropy.io import fits
from numpy import array, median, zeros, where, append, empty
from sys import argv, exit
from os import path, remove, listdir
import fnmatch
import glob
#import logutils
import statistics #Street's statistics functions
#import fitsutilsA
import warnings

#####

#####
# FUNCTIONS
def create_image_list_from_path(camera, date, caltype, startorend):

    list_of_images = listdir("{0}/{1}/{2}/{3}".format(camera,date,caltype,startorend))
    number_of_images = len(list_of_images)

    return list_of_images, number_of_images


def create_zeroed_image(camera, date, caltype, startorend, imagelist, number_of_images):

    if fnmatch.fnmatch(imagelist[0], "*00.fits"): # b00, d00, f00 Get dimensions from first image.
        init_image = fits.open("{0}/{1}/{2}/{3}/{4}".format(camera, date, caltype, startorend, imagelist[0]))
        Y_axis = init_image[0].header['NAXIS2']
        X_axis = init_image[0].header['NAXIS1']
        init_image.close()

        array_of_zeros = zeros((number_of_images,Y_axis,X_axis))

    return array_of_zeros


def assemble_image_cube(camera, date, caltype, startorend, imagelist, image_array):
    
    for image_index, image_name in enumerate(imagelist):

        if fnmatch.fnmatch(image_name, "*d00.fits"): # b00, d00, f00
            image = fits.open("{0}/{1}/{2}/{3}/{4}".format(camera, date, caltype, startorend, image_name))   # open FITS file
            imageheader = image[0].header
            imagedata = image[0].data
            image.close()

            image_array[image_index] = imagedata

    return imageheader, image_array


def outputFITS(imagedata,headerdict,filename): # copied from Street's SinistroToolKit
    '''Output data array to a FITS frame'''

    newhdu = fits.PrimaryHDU(imagedata)
    newhdulist = fits.HDUList([newhdu])
    for key,value in headerdict.items():
        newhdulist[0].header.update(key,value)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        newhdulist.writeto(filename,clobber=True,output_verify='ignore')
    newhdulist.close()

    return 0

def makebpm(MasterFile,nsigma_hi,nsigma_lo):
    '''Function to make a Bad Pixel Mask from a master flatfield'''

    # Read in the master flat field:
    if path.isfile(MasterFile) == False:
        print 'ERROR: Cannot find input FITS file ',MasterFile
        exit()
    image = fits.open(MasterFile)
    imageheader = image[0].header
    imagedata = image[0].data
    image.close()

    # Use smaller region to calculate mean
    stats_data = imagedata[5:1020,5:1530]

    # Compute statistics on the pixel values of the frame:
    (image_mean,image_sigma) = statistics.calcRMSclip2D(stats_data,1.0,0) # dataset, sigClip, iters

    print('Image mean and stddev: {0}, {1:.3f}'.format(image_mean,image_sigma))
    #print('Pix value at 1568, 1056: {0}'.format(imagedata[1055,1567]))

    # Identify all pixels above and below nsigma:
    thresh_hi = image_mean+float(nsigma_hi)*image_sigma
    thresh_lo = image_mean-float(nsigma_lo)*image_sigma
    print 'Applying thresholds HI=',thresh_hi,'ADU, LO=',thresh_lo

    idx1 = where(imagedata > thresh_hi)
    print("idx1 = {0}".format(idx1))
    #print("len(idx1) = {0}".format(len(idx1)))
    print("len(idx1[0]) = {0}".format(len(idx1[0])))
    print("len(idx1[1]) = {0}".format(len(idx1[1])))

    idx2 = where(imagedata < thresh_lo)

    # Generate FITS image of the Bad Pixel Mask:
    bpm = zeros([imageheader['NAXIS2'],imageheader['NAXIS1']])
    bpm[idx1] = 1.0
    bpm[idx2] = 1.0

    return bpm 


#####
# COMMAND LINE
if __name__ == '__main__':

    if len(argv) < 7:
        print 'Call sequence: '
        print 'python makebpm.py [Master Flat Field] [nsigma_hi] [nsigma_lo]'
        print 'where nsigma_hi,lo = float = factor of the image sigma to use as hi and lo pixel thresholds'
        exit()
    else:
        camera = argv[1]
        date = argv[2]
        caltype = argv[3]
        startorend = argv[4]
        nsigma_hi = float(argv[5])
        nsigma_lo = float(argv[6])

        imagelist, imagenumber = create_image_list_from_path(camera, date, caltype, startorend)

        zeroed_image = create_zeroed_image(camera, date, caltype, startorend, imagelist, imagenumber)

        imageheader, image_array = assemble_image_cube(camera, date, caltype, startorend, imagelist, zeroed_image)

        #median_image = median(image_array, axis=0)

        iexec = outputFITS(median(image_array, axis=0),imageheader,'{0}_med.fits'.format(camera))

        median_file = ('{0}_med.fits'.format(camera))
        
        bpm = makebpm(median_file,nsigma_hi,nsigma_lo)

        #print 'Output BPM ',median_file.replace('.fits','_bpm.fits')

        iexec = outputFITS(bpm,imageheader,median_file.replace('.fits','_bpm.fits'))

