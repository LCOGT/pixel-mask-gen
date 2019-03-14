import astropy.io.fits as fits
import pytest
import bpm.image_utils as image_utils

def test_apply_header_value_to_all_extensions():
    hdu_list = fits.HDUList(hdus=[fits.PrimaryHDU(header=fits.Header({'EXPTIME':10})),
                                  fits.ImageHDU(),
                                  fits.ImageHDU()])

    image_utils.apply_header_value_to_all_extensions([hdu_list], 'EXPTIME')

    assert hdu_list[1].header['EXPTIME'] == 10
    assert hdu_list[2].header['EXPTIME'] == 10
