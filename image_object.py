class ImageObject:
    """
    A class representing a FITS image object.

    Put more information here

    """

    def __init__(self, image_data, image_headers):
        """

        :param image_data:
        :param image_headers:

        """
        self.image_data = image_data
        self.image_headers = image_headers


    def get_image_size(self):
        """Returns an (x,y) tuple where x represents the number of rows and y represents the number of columns.

        :return:

        """
        return self.image_data.shape

    def get_image_header(self, key=None):
        """

        :param key: The 'key'/title of the image header you woud like to retrieve. If None, then the entire header\
        dictionary will be represented

        :return: If Key=None, then dictionary object representing the headers of the fits image, or if a key is specified,\
        then the corresponding value for the header key will be returned.

        :rtype: dict or string

        """

        if key is not None:
            return self.image_headers[key]

        return self.image_headers


    def get_image_data(self):
        """
        Self explanatory for now
        :return:

        """
        return self.image_data
