# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

WORKDIR /app

RUN apt update
RUN apt-get -y install build-essential

RUN apt -y install git-all

RUN mkdir ${HOME}/.ssh
RUN ssh-keyscan -t rsa github.com >> ${HOME}/.ssh/known_hosts

RUN git clone https://github.com/seansolari/expam.git
WORKDIR /app/expam

RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r requirements.txt
RUN python3 setup.py install

CMD ["expam", "--help"]
