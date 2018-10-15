from setuptools import setup, find_packages

setup(name='automated-pixel-mask-generation',
      author=['Curtis McCully', 'Raleigh Littles', 'Matt Daily'],
      version='0.1.0',
      packages=find_packages(),
      setup_requires=['pytest-runner'],
      install_requires=['astropy', 'lcogt-logging'],
      tests_require=['pytest', 'pytest-cov'])
