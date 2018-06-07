###################################
# IMPORTED FUNCTIONS
#import pyfits
from astropy.io import fits
import warnings

###############################
# OUTPUTFITS
def outputFITS(imagedata,headerdict,filename):
    '''Output data array to a FITS frame'''
    
    newhdu = fits.PrimaryHDU(imagedata)
    newhdulist = fits.HDUList([newhdu])
    for key,value in headerdict.items():
      	#newhdulist[0].header.update(key,value)
      	newhdulist[0].header.set(key,value)
    with warnings.catch_warnings():
      	warnings.simplefilter('ignore')
      	newhdulist.writeto(filename,clobber=True,output_verify='ignore')
    newhdulist.close()
    
    return 0
