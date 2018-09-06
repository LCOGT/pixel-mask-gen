
import src.image_processing as image_processing
import numpy as np
import astropy.io.fits as fits
import fractions #will be going away
import logging

logger = logging.getLogger(__name__)

def test_apply_bias_processing():
    test_image = np.random.rand(100,100)
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0,10)])
    test_image[tuple(bad_pixels.T)] = 10

    bias_mask = image_processing.apply_bias_processing([fits.ImageHDU(data=test_image)])

    assert np.shape(bias_mask) == np.shape(test_image)
    assert np.array_equal(np.argwhere(bias_mask == 1), np.argwhere(test_image == 10))

def test_apply_darks_processing():
    hdr = fits.Header()
    hdr['BIASSEC'] = '[1:5, 1:5]'
    hdr['EXPTIME'] = '10'

    test_image = np.round(30 + np.random.rand(100,100))
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0,10)])
    test_image[tuple(bad_pixels.T)] = 1000

    dark_mask = image_processing.apply_darks_processing([fits.ImageHDU(data=test_image,
                                                                       header=hdr)])
    assert np.shape(dark_mask) == np.shape(test_image)
    assert np.array_equal(np.argwhere(dark_mask == 1), np.argwhere(test_image == 1000))

def test_apply_flats_processing():
    hdr = fits.Header()
    hdr['FILTER'] = 'w'

    base_flat_image = np.random.randn(100,100) + 22000
    base_flat_std = np.std(base_flat_image)
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0,10)])
    base_flat_image[tuple(bad_pixels.T)] = 22000 + 8*base_flat_std

    flat_list = [fits.ImageHDU(header=hdr, data=base_flat_image * (0.9 ** index)) for index in range(1,6)]
    flat_mask = image_processing.apply_flats_processing(flat_list)

    assert np.shape(flat_mask) == np.shape(base_flat_image)
    assert np.array_equal(np.argwhere(flat_mask == 1),
                          np.argwhere(base_flat_image == 22000 + 8*base_flat_std))


def test_extract_center_fraction_region():
    test_image = np.zeros((100,100))
    test_image[24:75, 24:75] = 1

    assert np.array_equal(image_processing.extract_center_fraction_region(test_image, fractions.Fraction(1,4)),
                          np.ones((50,50)))

def test_combine_image_masks():
    mask_1 = np.array([1,0,0,1]).reshape(2,2)
    mask_2 = np.array([0,1,1,0]).reshape(2,2)

    combined_mask = image_processing.combine_image_masks([mask_1, mask_2])

    assert np.array_equal(combined_mask, np.ones((2,2)))

def test_extract_coordinates_from_header_string():
    test_header_string_1 = '[3100:3135, 1:2048]'
    test_header_string_2 = '[3100:3135,1:2048]'

    assert image_processing.extract_coordinates_from_header_string(test_header_string_1) ==\
           [3100, 3135, 1, 2048]

    assert image_processing.extract_coordinates_from_header_string(test_header_string_2) ==\
           [3100, 3135, 1, 2048]

def test_apply_frequency_thresholding_on_masked_arrays():
    test_array = np.array([5,1,2,4]).reshape(2,2)

    assert np.array_equal(image_processing.apply_frequency_thresholding_on_masked_arrays([test_array], 3),
                          np.array([True,False,False,True]).reshape(2,2))
