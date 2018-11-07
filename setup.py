from setuptools import setup, find_packages

setup(name='lco-bpm-maker',
      author=['Curtis McCully', 'Raleigh Littles', 'Matt Daily'],
      version='0.1.0',
      packages=find_packages(),
      setup_requires=['pytest-runner'],
      install_requires=['astropy', 'lcogt-logging'],
      tests_require=['pytest', 'pytest-cov'],
      entry_points={'console_scripts': ['create_bad_pixel_mask=bpm.generate_bpm:generate_bpm']})
