
# Project README

A Python script to generate bad pixel masks for a set of calibration images.

Full documentation is available at: `documentation/_build/html/index.html`

## Installation

#### Run
Run:
```bash
$ python3 setup.py build
$ python3 setup.py install
```

to build and install the project.

##### Development
For development, you way to use the `development`
option install of `build`.

Also, you may want to set up a virtual environment:
```bash
$ python3 -m virtualenv -p python3.6 virtual_environment
```

will set up a virtual environment and place its directories in `virtual_environment`

Once the virtual environment is set up, you must activate it:
```bash
$ source virtual_environment/bin/activate
```


#### Documentation building
Building documentation locally:
```bash
$ cd documentation/
$ make clean && make html
```
Then visit: `http://localhost:63342/pixel-mask-gen/documentation/_build/html/index.html` to see the documentation.
