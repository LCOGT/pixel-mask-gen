import src.image_processing as image_processing
import numpy as np
import astropy.io.fits as fits
import fractions #will be going away

def test_apply_bias_processing():
    hdr = fits.Header()

    test_image = np.round(np.random.normal(1, 1, (100, 100)))
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0, 10)])
    test_image[tuple(bad_pixels.T)] = 100

    bias_mask = image_processing.apply_bias_processing([fits.ImageHDU(data=test_image)])

    assert np.shape(bias_mask) == np.shape(test_image)
    np.testing.assert_array_equal((bias_mask == 1), (test_image == 100))

def test_apply_darks_processing():
    hdr = fits.Header()
    hdr['BIASSEC'] = '[1:5, 1:5]'
    hdr['EXPTIME'] = '10'

    test_image = np.round(np.random.normal(30, 1, (100,100)))
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0, 10)])
    test_image[tuple(bad_pixels.T)] = 1000

    dark_mask = image_processing.apply_darks_processing([fits.ImageHDU(header=hdr, data=test_image)])

    assert np.shape(dark_mask) == np.shape(test_image)
    np.testing.assert_array_equal((dark_mask == 1), (test_image == 1000))

def test_apply_flats_processing():
    hdr = fits.Header()
    hdr['FILTER'] = 'w'

    base_flat_image = np.random.normal(22000, 100, (100, 100))
    base_flat_std = np.std(base_flat_image)
    bad_pixels = np.array([[np.random.randint(100), np.random.randint(100)] for index in range(0, 10)])
    base_flat_image[tuple(bad_pixels.T)] = 22000 + 8*base_flat_std

    flat_list = [fits.ImageHDU(header=hdr, data=base_flat_image - (1000 * index)) for index in range(0, 5)]
    flat_mask = image_processing.apply_flats_processing(flat_list)

    assert np.shape(flat_mask) == np.shape(base_flat_image)
    np.testing.assert_array_equal((flat_mask == 1),
                                  (base_flat_image == 22000 + 8*base_flat_std))


def test_extract_center_fraction_region():
    test_image = np.zeros((100, 100))
    test_image[24:75, 24:75] = 1

    np.testing.assert_array_equal(image_processing.extract_center_fraction_region(test_image, fractions.Fraction(1, 4)),
                                  np.ones((50, 50)))

def test_combine_image_masks():
    mask_1 = np.array([1, 0, 0, 1]).reshape(2, 2)
    mask_2 = np.array([0, 1, 1, 0]).reshape(2, 2)

    combined_mask = image_processing.combine_image_masks([mask_1, mask_2])

    np.testing.assert_array_equal(combined_mask, np.ones((2, 2)))

def test_extract_coordinates_from_header_string():
    test_header_string_1 = '[3100:3135, 1:2048]'
    test_header_string_2 = '[3100:3135,1:2048]'

    assert image_processing.extract_coordinates_from_header_string(test_header_string_1) ==\
           [3100, 3135, 1, 2048]

    assert image_processing.extract_coordinates_from_header_string(test_header_string_2) ==\
           [3100, 3135, 1, 2048]

def test_apply_frequency_thresholding_on_masked_arrays():
    test_array = np.array([5, 1, 2, 4]).reshape(2, 2)

    np.testing.assert_array_equal(image_processing.apply_frequency_thresholding_on_masked_arrays([test_array], 3),
                          np.array([True, False, False, True]).reshape(2, 2))
