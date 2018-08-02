import io
import os
import sys

import logging
import lcogt_logging

from setuptools import find_packages, setup

"""
Taken from: https://github.com/kennethreitz/setup.py/blob/master/setup.py
"""

NAME = 'pixel-mask-gen'
DESCRIPTION = 'Bad pixel mask generator.'
URL = 'https://github.com/LCOGT/pixel-mask-gen'
EMAIL = 'rlittles@lco.global'
VERSION = None


AUTHOR = 'Raleigh Littles'
REQUIRES_PYTHON = '>=3.6'

REQUIRED = ['astropy','pyyaml', 'lcogt-logging']
EXTRAS = {

}

here = os.path.abspath(os.path.dirname(__file__))

about = {}
if not VERSION:
    with open(os.path.join(here, NAME, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION

# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    #long_description=long_description,
    #long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=('tests',)),
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],

    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ],
    # $ setup.py publish support.
)


def setup_custom_logger(name='pixel-mask-gen'):
    # Taken from: https://github.com/LCOGT/lcogt_logging/blob/master/example.py

    logger = logging.getLogger(name)
    helpful_info = {'example': True, 'pid': os.getpid()}
    formatter = lcogt-logging.LCOGTFormatter(extra_tags=helpful_info)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger

global LOGGER
LOGGER = setup_custom_logger()
