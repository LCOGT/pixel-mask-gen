import unittest
import numpy
import astropy.stats
import script
import random
import shutil
import logging

class TestFitsUtilities(unittest.TestCase):

    def test_writing_empty_array_as_fits(self):
        empty_array = numpy.empty([0, 0])
        with self.assertRaises(ValueError):
            script.output_to_FITS(empty_array, {}, 'filename.fits')

class TestsProcessingUtilities(unittest.TestCase):

    def test_median_filtering(self):
        random_array = numpy.random.rand(5,5)
        sigma_thresh = 1
        filtered_array = astropy.stats.sigma_clip(data=random_array, sigma=sigma_thresh)
        median = numpy.median(random_array)
        # The following array should be empty (it should have been filtered)
        bad_values = filtered_array[filtered_array > median + sigma_thresh]
        self.assertEqual(bad_values.size, 0)

    def test_adjacent_bad_pixels_detected(self):
        # construct a list of tuples where each tuple is an x,y coordinate
        # then at the end, take an existing coordinate and add 1 to the x-coordinate
        sim_bad_pixels_list = []
        for num in range(0,19):
            sim_bad_pixels_list.append((random.randint(1,100), random.randint(1,100)))

        existing_x_coord, existing_y_coord = sim_bad_pixels_list[random.randint(1,10)]

        sim_bad_pixels_list.append((existing_x_coord +1, existing_y_coord))

        # once neighboring bad pixel is added
        adjacent_bad_pixels = script.test_adjacent_pixels(sim_bad_pixels_list)
        self.assertGreaterEqual(adjacent_bad_pixels, 1)

class TestDirectorySettings(unittest.TestCase):

    def test_empty_image_list(self):
        image_list = []
        prefixes_array = ['b00', 'd00', 'f00']
        with self.assertRaises(ValueError):
            script.extract_data_from_files(image_list, prefixes_array)

    def test_invalid_config_file(self):
        with self.assertRaises(OSError):
            script.retrieve_image_directory_information('nonexistent_file.yml')


if __name__ == '__main__':
    unittest.main()



