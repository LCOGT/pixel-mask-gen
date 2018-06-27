# External imports
import unittest
import numpy
import os
import datetime
import subprocess
import re
import yaml
import sys
import random
import collections
import glob
import pdb
import astropy.io
import fractions

# Internal imports
import script
import image_object
import image_processing

class TestFitsFileIOUtilities(unittest.TestCase):

    def test_writing_empty_array_to_fits_raises_exception(self):
        empty_array = numpy.empty([0, 0])
        with self.assertRaises(ValueError):
            script.output_to_FITS(empty_array, {}, 'filename.fits')

    def test_empty_image_file_raises_exception(self):
    # make a list that contains the filename of one image
    # create an empty file with that filename
    # try to read in that filename, expect exception
        testing_prefixes = []
        empty_image_filename = os.path.join('test', 'empty_file.fits')

        subprocess.run(['touch', empty_image_filename])
        with self.assertRaises(IOError):
            script.extract_data_from_files([empty_image_filename], [''])

        subprocess.run(['rm', empty_image_filename])

class TestsImageProcessingUtilities(unittest.TestCase):

    def test_sigma_clipping(self):
        random_array = numpy.random.randint(5, size=(5,5))
        med, std = numpy.median(random_array), numpy.std(random_array)
        filtered_array, masked_indices, _ = image_processing.sigma_clip_individual(random_array, \
                                                                                   sigma_hi=1, \
                                                                                   sigma_low = 1)

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
        adjacent_bad_pixels = image_processing.test_adjacent_pixels(sim_bad_pixels_list)
        self.assertGreaterEqual(adjacent_bad_pixels, 1)


    def test_generate_mask_from_bad_pixels(self):
        # suppose there's bad pixels at the diagonals
        bad_pixel_locations = [(1,1), (2,2), (3,3)]
        bad_pixel_locations = [(1,1), (2,2), (3,3)]
        zero_array = numpy.zeros((4,4),dtype=bool)

        masked_array = script.generate_mask_from_bad_pixels(bad_pixel_locations, 4, 4)
        # masked array should have 3 bad pixels in it
        self.assertEqual(3, masked_array.sum())


    @unittest.skip('Not implemented')
    def test_bias_processing(self):
        # Do nearly the same process as flats
        test_image_num = 5
        test_image_size = 16,16

        test_images_list = []

        bad_pixel = random.randint(0, min(test_image_size)), \
                    random.randint(0, min(test_image_size))

        for image_index in range(test_image_num):
            random_array = numpy.random.randint(2, size=test_image_size)

            if random.getrandbits(1):
                random_array[bad_pixel] = sys.maxsize - 1

            test_images_list.append(image_object.ImageObject(random_array, {}))

        masked_indices = image_processing.flats_processing(test_images_list)

        self.assertEqual(masked_indices.size, 1)

    @unittest.skip('j')
    def test_flat_processing(self):
        # TODO: Use the setup/teardown methods for the processing types
        # Do a near identical process as flats
        test_image_num = 5
        test_image_size = 16,16

        test_images_list = []

        bad_pixel = random.randint(0, min(test_image_size) -1 ), \
                    random.randint(0, min(test_image_size) -1 )

        for image_index in range(test_image_num):
            random_array = numpy.random.randint(2, size=test_image_size)

            if random.getrandbits(1):
                random_array[bad_pixel] = sys.maxsize - 1

            test_images_list.append(image_object.ImageObject(random_array, {}))

        masked_indices = image_processing.flats_processing(test_images_list)

        self.assertEqual(masked_indices.size, 1)

    @unittest.skip('f')
    def test_dark_processing(self):
        # Create 5 image objects to represent 3 randomly generated images
        # Pick one pixel from one of the arrays that will be set to a known-bad value
        # Ensure it gets filtered

        test_image_num = 5
        test_image_size = 16,16
        # a list of numpy arrays that will represent our images?
        test_images_list = []

        # Pick a random pixel to use for testing
        bad_pixel = random.randint(0, min(test_image_size)), \
                    random.randint(0, min(test_image_size))

        for image_index in range(test_image_num):
            random_array = numpy.random.randint(2, size=test_image_size)
            image_header={'EXPTIME': random.uniform(1,10)}

            # Basically flip a coin and see if the bad pixel should be set
            if random.getrandbits(1):
                random_array[bad_pixel] = sys.maxsize - 1

            test_images_list.append(image_object.ImageObject(random_array, image_header))

        masked_indexes = image_processing.darks_processing(test_images_list)

        self.assertEqual(masked_indexes.size, 1)


    def test_center_region_extraction(self):
        """

        For the testing array of:

        1,  2,  3,  4
        5,  6,  7,  8
        9,  10, 11, 12
        13, 14, 15, 16

        The center quarter is 6, 7, 10, 11

        :return:
        """

        test_array = numpy.array([[1,2,3,4],
                                  [5,6,7,8],
                                  [9,10,11,12],
                                  [13,14,15,16]])

        expected_array = numpy.array([[6,7],
                                     [10,11]])

        extracted_array = image_processing.extract_center_fraction_region(test_array, fractions.Fraction(1,4))

        self.assertTrue(numpy.array_equal(expected_array, extracted_array))

class TestDateParsingAndPrefixes(unittest.TestCase):

    def test_empty_image_list_and_bad_prefix_arrays_raises_error(self):
        image_list = []
        prefixes_array = ['', '', '']
        with self.assertRaises(ValueError):
            script.extract_data_from_files(image_list, prefixes_array)

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


class TestConfigurationFileIssues(unittest.TestCase):

    def setUp(self):
        fake_config_filename = os.path.join('test', 'config.yml')
        self.fake_config_filename = fake_config_filename
        # check=True ensures an exception is thrown if the command returns a non-zero exit code,
        # which will prevent our tests in this class from running
        subprocess.run(['touch', fake_config_filename], check=True)

        with open('config.yml') as original_yml_file:
            configData = yaml.safe_load(original_yml_file)
            original_yml_file.close()
            self.original_configuration_data = configData

    def tearDown(self):
        subprocess.run(['rm', self.fake_config_filename], check=True)

    def test_empty_config_file_raises_error(self):
        # Create a fake empty yml file, then delete it
        # the error message has to say something like 'empty'
        with self.assertRaisesRegex(expected_exception=ValueError, expected_regex=re.compile(r'empty')):
            script.main(self.fake_config_filename)


    def test_invalid_yaml_structure_raises_error(self):
        # Create a yaml file that just contains one blank dictionary
        fake_data = {}
        with open(self.fake_config_filename, 'w') as yaml_file:
            yaml.dump(fake_data, yaml_file, default_flow_style=True)
            yaml_file.close()
        # since this is an error in parsing the yaml file, the error message should say something about this
        with self.assertRaisesRegex(expected_exception=ValueError, expected_regex=re.compile(r'parse')):
            script.main(self.fake_config_filename)

    def test_invalid_top_level_folder_halts_execution(self):
        # make a valid configuration file but with a bad top_level path

        # to avoid hardcoding in the configuration file structure, read the original configuration file in as a dict
        # then change the top directory parameter to something nonsense

        bad_config_data = self.original_configuration_data

        bad_config_data['directories']['top_directory'] = 'n/a'

        with open(self.fake_config_filename, 'w') as testing_yml_file:
            yaml.dump(bad_config_data, testing_yml_file, default_flow_style=True)
            testing_yml_file.close()

        with self.assertRaises(SystemExit):
            script.main(self.fake_config_filename)

    def test_invalid_image_folder_raises_exception(self):
        # make a valid configuration file with a real top_level_path thats empty
        # make a valid test configuration file, then make a real top_level path thats empty

        # to avoid hardcoding in the configuration file structure, read the original configuration file in as a dict
        # change the top_level directory (so that there isnt a folder conflict), then create the new directory

        fake_top_level_directory = os.path.join('test','top_level_directory')
        bad_config_data = self.original_configuration_data
        bad_config_data['directories']['top_directory']  = fake_top_level_directory

        with open(self.fake_config_filename, 'w') as testing_yml_file:
            yaml.dump(self.original_configuration_data, testing_yml_file, default_flow_style=True)
            testing_yml_file.close()

        subprocess.run(['mkdir', fake_top_level_directory], check=True)

        with self.assertRaises(FileNotFoundError):
            # the top level directory exists, but the image folder doesnt
            script.retrieve_image_directory_information(self.fake_config_filename, '69')

        subprocess.run(['rm', '-r', fake_top_level_directory], check=True)

    def test_valid_image_folder_with_no_images_raises_exception(self):
        # make a valid configuration file with a real top_level_path
        # in the top level path, make a folder for the date in the config file..
        # and then another folder for the camera prefix, but dont put anything in the folder

        # test that a ValueError is raised
        # delete the two directories you made, and the fake configuration file
        fake_top_level_directory = os.path.join('test','top_level_directory')
        bad_config_data = self.original_configuration_data
        bad_config_data['directories']['top_directory'] = fake_top_level_directory

        with open(self.fake_config_filename, 'w') as testing_yml_file:
            yaml.dump(bad_config_data, testing_yml_file, default_flow_style=True)
            testing_yml_file.close()

        # This tells you what folder to make
        date_string = script.parse_config_file_date(bad_config_data['directories'])

        fake_camera_id = '69'
        fake_image_path = os.path.join(fake_top_level_directory,bad_config_data['directories']['camera_prefix'] + fake_camera_id, date_string)
        subprocess.run(['mkdir', '-p', fake_image_path], check=True) # must use P flag to create nested directory

        with self.assertRaises(ValueError):
            script.retrieve_image_directory_information(self.fake_config_filename, '69')

        subprocess.run(['rm', '-r', fake_top_level_directory], check=True)


    def test_invalid_config_file_raises_error(self):
        with self.assertRaises(FileNotFoundError):
            script.retrieve_image_directory_information('nonexistent_file.yml', '69')

@unittest.skip('Testing')
class TestFullEndtoEnd(unittest.TestCase):

    def setUp(self):
        """To set up an integration test:
        1) Copy a valid yaml configuration file into the testing directory, but change the top_level directory into
        2) something that is within the test directory
        3) Create several random numpy arrays (to represent images) of fixed size, then select 1 arbitrary index (to rep-
        resent pixels) from each array that will be set to INTEGER_MAX, to represent a stuck pixel.
         4) Write these arrays to fits files and store the fits files in the path where the program will be expecting them
         (see earlier)

         The actual test will just consist of testing that the pixels you set to INT_MAX earlier were actually removed
        """
        # Set up an empty configuration yaml file in the test directory
        TestConfigurationFileIssues.setUp(self)

        testing_config_data = self.original_configuration_data
        self.fake_top_level_directory = os.path.join('test', 'top_level_directory')
        testing_config_data['directories']['top_directory'] = self.fake_top_level_directory
        testing_config_data['statistics']['threshold'] = 95
        # the prefixes list needs to change too
        new_prefix_list = ['b69', 'd69', 'f69']
        #  its okay to hard code these in because these are the only calibration types that will ever exist
        old_prefix_words = ['bias', 'dark', 'flats']
        for index, value in enumerate(new_prefix_list):
            testing_config_data['directories'][str(old_prefix_words[index]) + '_prefix'] = value

        print('Initiating YAML file.')

        with open(self.fake_config_filename, 'w') as testing_yml_file:
            yaml.dump(testing_config_data, testing_yml_file, default_flow_style=True)
            testing_yml_file.close()

        print('Test YAML file initated.')

        date_string = script.parse_config_file_date(testing_config_data['directories'])

        fake_camera_id = '69'
        fake_image_path = os.path.join(self.fake_top_level_directory, testing_config_data['directories']['camera_prefix'] + fake_camera_id, date_string)

        subprocess.run(['mkdir', '-p', fake_image_path], check=True)

        # Once the path is setup, you know where to put the fits files you will create, so get ready to start creating them
        test_images = 12
        self.bad_pixel_locations = []
        self.image_dimensions = 2000, 2000

        max_pixel_val = 2
        print("Preparing to generate fake images into directory: {0}".format(fake_image_path))
        # select a random tuple to serve as the coordinates 'marked' pixel
        # use min to prevent an indexerror
        bad_pixel_index = ( random.randint(0, min(self.image_dimensions[0], self.image_dimensions[1])), \
                            random.randint(0, min(self.image_dimensions[0], self.image_dimensions[1])) )

        print("The fake bad pixel is in location: {0}".format(bad_pixel_index))

        for image_number in range(0, test_images):
            image_data = numpy.random.randint(low=0, high=max_pixel_val, size=(self.image_dimensions[0], self.image_dimensions[1]))

            bad_pixel_value = sys.maxsize - 1
            self.bad_pixel_locations.append(tuple(bad_pixel_index))
            image_data[bad_pixel_index] = bad_pixel_value
            current_prefix_num = image_number % len(new_prefix_list)
            test_image_filename = os.path.join(fake_image_path,  str(image_number) + "-random-test-{0}.fits".format(new_prefix_list[current_prefix_num]))

            # Write the image_data array to a .fits file
            new_hdu = astropy.io.fits.PrimaryHDU(image_data.astype(numpy.uint8))
            new_hdu_list = astropy.io.fits.HDUList([new_hdu])

            new_hdu_list[0].header.set('OBSTYPE', 'BPM')
            new_hdu_list.writeto(test_image_filename,overwrite=False,output_verify='exception')
            new_hdu_list.close()

        print('Completed setup for integration test.')
        # Create a FITS file that looks something like: path/to/images/here/<XY>_bpm-test.fits' that will be the sample images

    def tearDown(self):
        # Delete the testing configuration yaml file
        # recursively delete the fake testing top directory
        # Delete the debug directory?
        subprocess.run(['rm', self.fake_config_filename], check=True)
        subprocess.run(['rm', '-rf', self.fake_top_level_directory])


    def test_full_integration_test(self):
        # Now that all the sample images are created, call the main script?
        script.main(self.fake_config_filename)
        # find if any bad pixel appeared with a frequency that exceeds the threshold amount (for each of the three thirds)

        # Even though the threshold value is configurable, we hardcoded it to 95%
        thresholded_array = [key for (key, value) in collections.Counter(self.bad_pixel_locations).items() \
                             if (value >= 0.95 * len(self.bad_pixel_locations))]

        # once you have the final array, generate the mask
        masked_array = numpy.zeros(self.image_dimensions, dtype=bool)
        for coordinates in self.bad_pixel_locations:
            masked_array[coordinates] = True

        # after the main script has run, read out the bpm file produced by that mask and test it against the saved mask earlier, they should match?
        # the final bpm is a .fits file in the main directory, so search for it?
        fits_files_in_dir = glob.glob('*.fits')
        if len(fits_files_in_dir) < 2:
            combined_bpm_mask_filename = os.path.abspath(fits_files_in_dir[0])

        combined_bpm_data, combined_bpm_headers, combined_bpm_size = script.read_individual_fits_file(combined_bpm_mask_filename)

        # The mask here should have the same amount of nonzero elements as your mask did earlier, to ensure that the
        # correct amount of bad pixels were detected

        self.assertEqual(numpy.count_nonzero(combined_bpm_data), numpy.count_nonzero(masked_array))


if __name__ == '__main__':
    unittest.main()



