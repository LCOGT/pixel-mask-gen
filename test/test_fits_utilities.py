import os
import numpy
import unittest
import subprocess
import re

import fits_utilities


class TestFitsFileIOUtilities(unittest.TestCase):

    def test_writing_empty_array_to_fits_raises_exception(self):
        empty_array = numpy.empty([0, 0])
        with self.assertRaisesRegex(expected_exception=ValueError, expected_regex=re.compile(r'empty')):
            fits_utilities.output_to_FITS(empty_array, {}, 'filename.fits')

    def test_empty_image_file_raises_exception(self):
        # make a list that contains the filename of one image
        # create an empty file with that filename
        # try to read in that filename, expect exception
        empty_image_filename = os.path.join('test', 'empty_file.fits')

        subprocess.run(['touch', empty_image_filename])
        with self.assertRaisesRegex(expected_exception=ValueError, expected_regex=re.compile(r'blank')):
            fits_utilities.read_individual_fits_file('')

        subprocess.run(['rm', empty_image_filename])


if __name__ == '__main__':
    unittest.main()
