import src.image_processing as image_processing
import numpy as np
import astropy.io.fits as fits
import astropy.stats

def test_process_bias_frames():
    bad_pixel_locations = generate_bad_pixel_locations(94, 100, 10)
    bias_frames = [generate_test_bias_frame(bad_pixel_locations) for index in range(0,10)]

    bias_mask = image_processing.process_bias_frames(bias_frames, mad_threshold=8)
    flagged_pixels = np.where(bias_mask == True)

    assert np.shape(bias_mask) == np.shape(bias_frames[0].data)
    assert set(bad_pixel_locations[0]) == set(flagged_pixels[0])
    assert set(bad_pixel_locations[1]) == set(flagged_pixels[1])

def test_process_dark_frames():
    bad_pixel_locations = generate_bad_pixel_locations(94, 100, 10)
    dark_frames = [generate_test_dark_frame(bad_pixel_locations) for index in range(0,10)]

    dark_mask = image_processing.process_dark_frames(dark_frames)
    flagged_pixels = np.where(dark_mask==True)

    assert np.shape(dark_mask) == np.shape(dark_frames[0].data)
    assert set(flagged_pixels[0]) == set(bad_pixel_locations[0])
    assert set(flagged_pixels[1]) == set(bad_pixel_locations[1])

def test_process_flat_frames():
    hdr = fits.Header()
    hdr['FILTER'] = 'w'
    hdr['BIASSEC'] = '[95:100, 1:100]'
    hdr['TRIMSEC'] = '[1:94, 1:100]'

    base_image_mean = 22000
    base_image_std = 1000
    bad_pixel_locations = generate_bad_pixel_locations(94, 100, 10)

    flat_frames = [generate_test_flat_frame(bad_pixel_locations,
                                            base_image_mean - 2000*index,
                                            base_image_std - 50*index) for index in range(0,10)]

    flat_mask = image_processing.process_flat_frames(flat_frames, mad_threshold=8)
    flagged_pixels = np.where(flat_mask == True)

    assert np.shape(flat_mask) == np.shape(flat_frames[0].data)
    assert set(bad_pixel_locations[0]) == set(flagged_pixels[0])
    assert set(bad_pixel_locations[1]) == set(flagged_pixels[1])


def test_get_slices_from_image_section():
    test_header_string_1 = '[3100:3135, 1:2048]'
    test_header_string_2 = '[3100:3135,1:2048]'

    assert image_processing.get_slices_from_image_section(test_header_string_1) ==\
           (slice(0, 2048, 1), slice(3099, 3135, 1))

    assert image_processing.get_slices_from_image_section(test_header_string_2) ==\
           (slice(0, 2048, 1), slice(3099, 3135, 1))


def generate_test_bias_frame(bad_pixel_locations):
    hdr = fits.Header()
    hdr['BIASSEC'] = '[95:100, 1:100]'

    overscan_slices = image_processing.get_slices_from_image_section(hdr['BIASSEC'])

    bias_frame = np.round(np.random.normal(1000, 50, (100, 100)))
    bias_frame[bad_pixel_locations] = np.random.normal(5000, 3000)
    bias_frame[overscan_slices] = np.round(np.random.normal(1000, 50))

    return fits.ImageHDU(data=bias_frame, header=hdr)

def generate_test_flat_frame(bad_pixel_locations, image_mean, image_std):
    hdr = fits.Header()
    hdr['FILTER'] = 'w'
    hdr['BIASSEC'] = '[95:100, 1:100]'
    hdr['TRIMSEC'] = '[1:94, 1:100]'

    overscan_slices = image_processing.get_slices_from_image_section(hdr['BIASSEC'])

    flat_frame = np.round(np.random.normal(image_mean, image_std, (100,100)))
    flat_frame[bad_pixel_locations] = np.round(np.random.normal(image_mean, 20*image_std))
    flat_frame[overscan_slices] = np.round(np.random.normal(1000, 50))

    return fits.ImageHDU(data=flat_frame, header=hdr)

def generate_test_dark_frame(bad_pixel_locations):
    hdr = fits.Header()
    hdr['BIASSEC'] = '[95:100, 1:100]'
    hdr['EXPTIME'] = '10.0'

    overscan_region_coordinates = image_processing.get_slices_from_image_section(hdr['BIASSEC'])

    dark_frame = np.round(np.random.normal(30, 5, (100,100)))
    dark_frame[bad_pixel_locations] = 500

    return fits.ImageHDU(data=dark_frame, header=hdr)

def generate_bad_pixel_locations(x_limit, y_limit, num_pixels):
    bad_pixels_y = np.array([np.random.randint(y_limit) for index in range(0, num_pixels)])
    bad_pixels_x = np.array([np.random.randint(x_limit) for index in range(0, num_pixels)])
    bad_pixels = tuple((bad_pixels_y, bad_pixels_x))

    return bad_pixels
