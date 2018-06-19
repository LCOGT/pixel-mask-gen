FROM python:3

WORKDIR /app/

MAINTAINER Raleigh Littles <rlittles@lco.global>

LABEL "repo-name"="pixel-mask-gen"

LABEL "repo-link"="https://github.com/LCOGT/pixel-mask-gen"

ADD script.py /

RUN pip install pipenv

RUN pipenv --three

RUN pipenv install -r app/requirements.txt

CMD ["pipenv", "run", "python3", "app/script.py", "app/config.yml"]
