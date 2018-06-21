
# Project README

A Python script to generate bad pixel masks for a set of calibration images.

Full documentation is available at: `documentation/_build/html/index.html`

## Pre-requisites

* Python 3.5

## Installation

### Docker

Run:

```bash
$ sudo docker build --no-cache -t pixel-mask-gen .

$ sudo docker run pixel-mask-gen

```

### Local

Instructions on how to set up virtual environment, install packages into virtualenvironment, run the tests, and run the code.
#### Run
Run:
```bash
$ python3 -m virtualenv -p python3.5 .
$ source bin/activate
$ pip install -r requirements.txt
$ python3 -m unittest -v test/test_script.py
$ python3 script.py config.yml
```

#### Testing
Getting test coverage:
```bash
$  coverage run -m unittest discover test/ --verbose --locals
$ coverage html

```
Then visit `/htmlcov/index.html` to see the test coverage.
#### Documentation
Building documentation:
```bash
$ cd documentation/ && make html && cd ..
```
Then visit: `http://localhost:63342/pixel-mask-gen/documentation/_build/html/index.html` to see the documentation.
