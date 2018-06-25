FROM python:3

MAINTAINER Raleigh Littles <rlittles@lco.global>

LABEL "repo-name"="pixel-mask-gen"

LABEL "repo-link"="https://github.com/LCOGT/pixel-mask-gen"

# you have to ADD every file that your code will need to run (?)
ADD script.py /

ADD requirements.txt /

ADD config.yml /

ADD test/test_script.py test/test_script.py

ADD debug/ debug/

ADD documentation/ documentation/

# Install requirements
RUN ["pip" ,"install", "-r", "requirements.txt"]

# Run tests
RUN ["python3", "-m", "unittest", "-v", "test/test_script.py"]


# CMD tells docker what to do when you actually do "docker run <my-image>".
# THERE CAN ONLY BE ONE OF THESE PER DOCKERFILE
CMD ["python3", "script.py", "config.yml"]
