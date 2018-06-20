import unittest
import numpy
import astropy.stats
import script
import random
import datetime
import subprocess
import re
import yaml

class TestFitsUtilities(unittest.TestCase):

    def test_writing_empty_array_as_fits(self):
        empty_array = numpy.empty([0, 0])
        with self.assertRaises(ValueError):
            script.output_to_FITS(empty_array, {}, 'filename.fits')


    def test_reading_empty_fits_file(self):
        # create an empty fits file
        # add its filename to a list
        # call extract_data_from_files() and expect ValueError

        # delete empty fits file

class TestsProcessingUtilities(unittest.TestCase):

    def test_median_filtering(self):
        random_array = numpy.random.randint(5, size=(5,5))
        med, std = numpy.median(random_array), numpy.std(random_array)
        print("med:{0}, std:{1}".format(med, std))
        filtered_array, masked_indices, _ = script.filter_individual(random_array, sigma_hi=1, sigma_low = 1)
        # we dont know how many iterations the algorithm will need to do, but we do know that its HIGHLY likely that at
        # least one entry out of the 25 will need to be filtered
        self.assertGreater(masked_indices.size, 1)

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


    def test_generate_mask_from_bad_pixels(self):
        # suppose there's bad pixels at the diagonals
        bad_pixel_locations = [(1,1), (2,2), (3,3)]
        zero_array = numpy.zeros((4,4),dtype=bool)

        masked_array = script.generate_mask_from_bad_pixels(bad_pixel_locations, 4, 4)
        # masked array should have 3 bad pixels in it
        self.assertEqual(3, masked_array.sum())

class TestDirectorySettings(unittest.TestCase):

    def test_empty_image_list_and_bad_prefix_arrays(self):
        image_list = []
        prefixes_array = ['', '', '']
        with self.assertRaises(ValueError):
            script.extract_data_from_files(image_list, prefixes_array)

    def test_invalid_config_file(self):
        with self.assertRaises(OSError):
            script.retrieve_image_directory_information('nonexistent_file.yml', '69')

    def test_empty_config_file(self):
        # Create a fake empty yml file, then delete it
        fake_filename = 'test/config.yml'
        subprocess.run(["touch", fake_filename])
        # the error message has to say something like 'empty'
        with self.assertRaisesRegex(expected_exception=ValueError, expected_regex=re.compile(r'empty')):
            script.main(fake_filename)

        subprocess.run(['rm', fake_filename])

    def test_invalid_yaml_structure_for_date_parsing(self):
        # Create a yaml file that just contains one blank dictionary
        fake_filename = 'test/config.yml'
        subprocess.run(["touch", fake_filename])

        fake_data = {}
        with open(fake_filename, 'w') as yaml_file:
            yaml.dump(fake_data, yaml_file, default_flow_style=True)
        # since this is an error in parsing the yaml file, the error message should say something about this
        with self.assertRaisesRegex(expected_exception=ValueError, expected_regex=re.compile(r'parse')):
            script.main(fake_filename)

        subprocess.run(['rm', fake_filename])

    # TODO: Remove code duplication in next two tests
    def test_bad_relative_date_returns_yesterdays_date(self):
        directory_info_dict = {'top_directory': '/notapplicable/',
                               'camera_prefix': 'xx', 'time': 'start',
                               'bias_prefix': 'b00', 'dark_prefix': 'd00', 'flats_prefix': 'f00'}

        directory_info_dict['date'] = {'exact': False, 'relative': '35th of nevuary'}

        date_string = script.parse_config_file_date(directory_info_dict)
        yesterdays_date = (datetime.datetime.today() - datetime.timedelta(days=1)).strftime('%Y%m%d')
        self.assertEqual(yesterdays_date, date_string)

    def test_bad_exact_date_returns_yesterdays_date(self):
        directory_info_dict = {'top_directory': '/notapplicable/',
                               'camera_prefix': 'xx', 'time': 'start',
                               'bias_prefix': 'b00', 'dark_prefix': 'd00', 'flats_prefix': 'f00'}

        directory_info_dict['date'] = {'exact': 'baddate', 'relative': '35th of nevuary'}

        yesterdays_date = (datetime.datetime.today() - datetime.timedelta(days=1)).strftime('%Y%m%d')
        date_string = script.parse_config_file_date(directory_info_dict)
        self.assertEqual(yesterdays_date, date_string)

    def test_invalid_image_folder(self):
        # make a valid configuration file but with a bad top_level path
        pass
        # then delete the test configuration file

    def valid_image_folder_with_no_images(self):
        # make a valid test configuration file, then make a real top_level path thats empty
        pass
        # then delete both the new path and the testing configuration file


if __name__ == '__main__':
    unittest.main()



