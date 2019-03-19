import astropy.io.fits as fits
import pytest
import numpy as np
import bpm.image_utils as image_utils

def test_apply_header_value_to_all_extensions():
    hdu_list = fits.HDUList(hdus=[fits.PrimaryHDU(header=fits.Header({'EXPTIME':10})),
                                  fits.ImageHDU(),
                                  fits.ImageHDU()])

    image_utils.apply_header_value_to_all_extensions([hdu_list], 'EXPTIME')

    assert hdu_list[1].header['EXPTIME'] == 10
    assert hdu_list[2].header['EXPTIME'] == 10

def test_get_sci_extensions_from_amplifier():
    hdu_list = fits.HDUList(hdus=[fits.PrimaryHDU(),
                                  fits.ImageHDU(data=np.array([0,0,0,0]),
                                                header=fits.Header({'EXTNAME': 'SCI',
                                                                    'EXTVER': 1})),
                                  fits.ImageHDU(data=np.array([1,1,1,1]),
                                                header=fits.Header({'EXTNAME': 'SCI',
                                                                    'EXTVER': 2}))])

    frames = image_utils.get_extensions_by_name(hdu_list, 'SCI')
    amplifier_1_frames = image_utils.get_sci_extensions_from_amplifier(frames, 1)
    amplifier_2_frames = image_utils.get_sci_extensions_from_amplifier(frames, 2)
    amplifier_3_frames = image_utils.get_sci_extensions_from_amplifier(frames, 3)

    assert len(amplifier_1_frames) == 1
    assert amplifier_1_frames[0].data.all() == np.array([0,0,0,0]).all()
    assert len(amplifier_2_frames) == 1
    assert amplifier_2_frames[0].data.all() == np.array([1,1,1,1]).all()
    assert amplifier_3_frames == []

def test_sort_frames_by_header_values():
    hdu_list = fits.HDUList(hdus=[fits.PrimaryHDU(),
                                  fits.ImageHDU(data=np.array([0,0,0,0]),
                                                header=fits.Header({'EXTNAME': 'SCI',
                                                                    'EXTVER': 1})),
                                  fits.ImageHDU(data=np.array([1,1,1,1]),
                                                header=fits.Header({'EXTNAME': 'SCI',
                                                                    'EXTVER': 2}))])

    frames = image_utils.get_extensions_by_name(hdu_list, 'SCI')
    frames_sorted = image_utils.sort_frames_by_header_values(frames, 'EXTVER')

    assert list(frames_sorted.keys()) == [1, 2]
    assert frames_sorted[1][0].data.all() == np.array([0,0,0,0]).all()
    assert frames_sorted[2][0].data.all() == np.array([1,1,1,1]).all()
