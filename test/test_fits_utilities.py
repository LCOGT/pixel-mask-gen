import unittest
import numpy

import src2.fits_utilities

class TestFitsUtilities(unittest.TestCase):

    def test_writing_to_fits(self):
        # Generate empty array
        empty_array = numpy.empty([0, 0])
        with self.assertRaises(BaseException):
            src2.fits_utilities.output_to_FITS(empty_array, {}, 'filename.fits')



if __name__ == '__main__':
    unittest.main()



