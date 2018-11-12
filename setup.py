from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='lco-bpm-maker',
      author='Curtis McCully, Raleigh Littles, Matt Daily',
      description='LCO bad pixel mask creator',
      long_description=long_description,
      version='0.1.0',
      packages=find_packages(),
      setup_requires=['pytest-runner'],
      install_requires=['numpy', 'astropy', 'lcogt-logging'],
      tests_require=['pytest', 'pytest-cov'],
      entry_points={'console_scripts': ['lco_bpm_maker=bpm.generate_bpm:generate_bpm']})
