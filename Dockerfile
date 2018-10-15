FROM python:3.6-slim-jessie

WORKDIR /lco/automated-bad-pixel-mask

COPY . /lco/automated-bad-pixel-mask

RUN python setup.py install

ENTRYPOINT ["python", "src/script.py"]
