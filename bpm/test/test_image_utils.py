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

def test_stack_frames_by_extver():
    frame_list_1 = fits.HDUList(hdus=[fits.ImageHDU(header=fits.Header({'EXTVER':1})),
                                      fits.ImageHDU(header=fits.Header({'EXTVER':2})),
                                      fits.ImageHDU(header=fits.Header({'EXTVER':3})),
                                      fits.ImageHDU(header=fits.Header({'EXTVER':4}))])

    frame_list_2 = fits.HDUList(hdus=[fits.ImageHDU(header=fits.Header({'EXTVER':1})),
                                      fits.ImageHDU(header=fits.Header({'EXTVER':2}))])

    sorted_frames = image_utils.stack_frames_by_extver([frame_list_1, frame_list_2])

    lengths = []
    extver_values = []
    for extver in sorted_frames.keys():
        lengths.append(len(sorted_frames[extver]))
        extver_values.extend([frame[0].header['EXTVER'] for frame in sorted_frames[extver]])

    assert lengths == [2, 2, 1, 1]
    assert extver_values == [1, 1, 2, 2, 3, 4]
    
