# pixel-mask-gen

A Python script to generate bad pixel masks for a set of calibration images.

Full documentation is available at: `documentation/_build/html/index.html`

# Pre-requisites

* Python 3.5

# Installation

## Docker

Run:

```bash
$ sudo docker build --no-cache -t pixel-mask-gen .

$ sudo docker run pixel-mask-gen

```

## Local

Run:
```bash
$ pip install -r requirements.txt
$ python3 -m unittest -v test/test_script.py
$ python3 script.py config.yml
```


