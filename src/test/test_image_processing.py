
import src.image_processing as image_processing
import numpy as np
import astropy.io.fits as fits

def test_apply_bias_processing():
    test_image = np.random.rand(100,100)
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0,10)])
    test_image[tuple(bad_pixels.T)] = 10

    bias_mask = image_processing.apply_bias_processing([fits.ImageHDU(data=test_image)])
    
    assert np.array_equal(np.argwhere(bias_mask == 1), np.argwhere(test_image == 10))

def test_apply_darks_processing():
    hdr = fits.Header()
    hdr['BIASSEC'] = '[1:5, 1:5]'
    hdr['EXPTIME'] = '10'

    test_image = 29.5 + np.random.rand(100,100)
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0,10)])
    test_image[tuple(bad_pixels.T)] = 1000

    dark_mask = image_processing.apply_darks_processing([fits.ImageHDU(data=test_image,
                                                                       header=hdr)])

    assert np.array_equal(np.argwhere(dark_mask == 1), np.argwhere(test_image == 1000))

def test_apply_flats_processing():
    pass

def test_extract_center_fraction_region():
    test_image = np.zeros((100,100))
    test_image[24:75, 24:75] = 1

    assert np.array_equal(image_processing.extract_center_fraction_region(test_image),
                          np.ones((50,50)))

def test_combine_image_masks():
    pass

def test_extract_coordinates_from_header_string():
    pass