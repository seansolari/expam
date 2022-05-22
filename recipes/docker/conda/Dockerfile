# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
RUN echo "conda activate expam" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# The code to run when container is started:
COPY entrypoint.sh ./
ENTRYPOINT ["./entrypoint.sh"]
