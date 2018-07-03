class ImageObject:
    """
    A class built for representing a FITS v4 image. Currently only supports single-extension FITS\
    images.

    """

    def __init__(self, image_data, image_headers):
        self.image_data = image_data
        self.image_headers = image_headers


    def get_image_size(self):
        """

         :return: a 2-element tuple where the first element represents the number of rows\
        and the second element represents the number of columns.

        """
        return self.image_data.shape

    def get_image_header(self, key=None):
        """

        :param key: The 'key'/title of the image header you woud like to retrieve. If None, then the entire header\
        dictionary will be returned.

        :return: If Key=None, then dictionary object representing the headers of the fits image, or if a key is specified,\
        then the corresponding value for the header key will be returned.

        :rtype: dict or string

        """

        if key is not None:
            return self.image_headers[key]

        return self.image_headers


    def get_image_data(self):
        """
        :return: the data stored in the first HDU as a numpy.array.

        For more information, see `the Astropy documentation <http://docs.astropy.org/en/stable/io/fits/>`_.

        """
        return self.image_data
