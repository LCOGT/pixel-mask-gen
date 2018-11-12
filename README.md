# pixel-mask-gen
Python utility to generate a bad pixel mask from a set of calibration images.

Authors: Matt Daily, Raleigh Littles, Curtis McCully

## Installation
### From PyPi
This package is available via PyPi, and can be installed via pip:

`pip3 install lco-bpm-maker`

### From Github
To install the tool, clone this repository and run:

```
cd pixel-mask-gen
python3 setup.py install
```

## Tests
To run the unit tests, simply run:

`python3 setup.py test`

## Usage
Once you've installed the tool, it can be run simply by:

`lco_bpm_maker`

```
usage: lco_bpm_maker [-h] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                     [--dark-current-threshold DARK_CURRENT_THRESHOLD]
                     [--flat-sigma-threshold FLAT_SIGMA_THRESHOLD]
                     [--bias-sigma-threshold BIAS_SIGMA_THRESHOLD]
                     input_directory output_directory

Create a bad pixel mask from a set of calibration frames.

positional arguments:
  input_directory       Input directory of calibration images

  output_directory      Output directory for bad pixel mask

optional arguments:
  -h, --help            show this help message and exit

  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level to be displayed

  --dark-current-threshold DARK_CURRENT_THRESHOLD
                        Threshold for pixel dark current when flagging bad
                        pixels in dark frames. Pixels above this will be
                        flagged. Default = 35 [electrons/second]

  --flat-sigma-threshold FLAT_SIGMA_THRESHOLD
                        Number of standard deviations from the median of the
                        combined flat image for a pixel to be flagged. Default = 10

  --bias-sigma-threshold BIAS_SIGMA_THRESHOLD
                        Number of standard deviations from the median of the
                        combined bias image for a pixel to be flagged. Default = 10
```
