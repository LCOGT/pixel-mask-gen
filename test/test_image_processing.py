import numpy
import random
import fractions
import unittest
import sys

import src.image_processing
import src.script
import src.image_object


class TestsImageProcessingUtilities(unittest.TestCase):

    def test_sigma_clipping(self):
        random_array = numpy.random.randint(5, size=(5,5))
        med, std = numpy.median(random_array), numpy.std(random_array)
        filtered_array, masked_indices, _ = src.image_processing.sigma_clip_individual(random_array, \
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
        adjacent_bad_pixels = src.image_processing.test_adjacent_pixels(sim_bad_pixels_list)
        self.assertGreaterEqual(adjacent_bad_pixels, 1)

    def test_generate_mask_from_bad_pixels(self):

        # Create a random array with bad pixels at the diagonal elements
        max_array_size = random.randint(5,15)

        zero_array = numpy.zeros((max_array_size, max_array_size),dtype=bool)

        bad_pixel_locations = [(i, i) for i in range(0, max_array_size)]

        for coords in bad_pixel_locations: zero_array[coords] = True

        masked_array = src.script.generate_mask_from_bad_pixels(bad_pixel_locations, max_array_size, max_array_size)
        # masked array should have 3 bad pixels in it
        self.assertEqual(len(bad_pixel_locations), masked_array.sum())

    def test_center_region_extraction(self):
        """
        For the testing array of:

        1,  2,  3,  4
        5,  6,  7,  8
        9,  10, 11, 12
        13, 14, 15, 16

        The center quarter is 6, 7, 10, 11
        """
        test_array = numpy.array([[1,2,3,4],
                                  [5,6,7,8],
                                  [9,10,11,12],
                                  [13,14,15,16]])

        expected_array = numpy.array([[6,7],
                                     [10,11]])

        extracted_array = src.image_processing.extract_center_fraction_region(test_array, fractions.Fraction(1,4))

        self.assertTrue(numpy.array_equal(expected_array, extracted_array))


class TestImageProcessingCoreFunctions(unittest.TestCase):

    def setUp(self):
        # Set up a set of random arrays that will serve as the test data
        num_of_test_images = 10

        # set a maximum image value for a pixel, 10 is a good number
        max_pixel_value = 100

        # to make the fractioning easy, use a factor of 4
        test_image_size = 16, 16

        self.test_images_list = []

        # randomly choose an index that will serve as the location of the bad pixel
        bad_pixel = ( random.randint(0, min(test_image_size) - 1), random.randint(0, min(test_image_size) - 1) )

        # This is only needed by the darks
        image_header={'EXPTIME': random.uniform(1,10)}

        self.expected_bad_pixel_count = 0

        for image_index in range(num_of_test_images):
            # this is the array that is going to represent an image
            random_array = numpy.random.randint(max_pixel_value, size=test_image_size)

            # 'flip' a coin to see if the current image is giong to have a bad pixel in it
            if random.getrandbits(1):
                # Pick the value that is going to actually be set in the bad pixel, dont use negative numbers since
                # our array will be converted to uint8
                # If the bad pixel is set, it needs to show up
                random_array[bad_pixel] = sys.maxsize - 1
                self.expected_bad_pixel_count += 1


            self.test_images_list.append(src.image_object.ImageObject(random_array, image_header))

        if len(self.test_images_list) != num_of_test_images:
            raise unittest.SkipTest("Test skipped, unable to create to testing images.")

    def test_bias_processing(self):
        masked_indices = src.image_processing.biases_processing(self.test_images_list, sigma_min=5, sigma_max=5)

        self.assertGreaterEqual(len(masked_indices), self.expected_bad_pixel_count)

    def test_flat_processing(self):
        masked_indices = src.image_processing.flats_processing(self.test_images_list)
        # flats processing wont contain duplicates in its bad pixel list, whereas our sample images
        # have duplicated pixels
        self.assertGreaterEqual(len(masked_indices), self.expected_bad_pixel_count)

    def test_dark_processing(self):
        masked_indexes = src.image_processing.darks_processing(self.test_images_list)

        self.assertGreaterEqual(len(masked_indexes), self.expected_bad_pixel_count)


if __name__ == '__main__':
    unittest.main()
