import unittest
import numpy
import astropy.stats
import script

class TestFitsUtilities(unittest.TestCase):

    def test_writing_empty_array_as_fits(self):
        empty_array = numpy.empty([0, 0])
        with self.assertRaises(BaseException):
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

    def test_rectangle_detected(self):
        pass

class TestDirectorySettings(unittest.TestCase):

    def test_empty_image_list(self):
        pass # should fail with exception

    def test_invalid_date(self):
        pass # should fail with exception

    def test_invalid_directory(self):
        pass # should fail with exception



if __name__ == '__main__':
    unittest.main()



