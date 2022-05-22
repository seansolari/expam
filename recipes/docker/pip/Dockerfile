# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

WORKDIR /app

# Install expam requirements.
RUN apt update
RUN apt-get -y install build-essential

# Install expam.
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install expam

RUN mkdir DATABASE
RUN mkdir INPUT
RUN mkdir OUTPUT

# Download the publicly available expam (test) database.
WORKDIR /app/DATABASE
RUN python3 -m pip install gdown
RUN gdown https://drive.google.com/uc?id=1KRrEvG5Sr28wvkEFW5CXKwcYqLHAu6f7

RUN tar -xvzf test.tar.gz
RUN rm test.tar.gz
RUN rm ._test

# Run database classification.
CMD [ "expam", "classify", "-db", "test", "-d", "/app/INPUT/reads", "-o", "/app/OUTPUT/my_run", "--taxonomy" ]
