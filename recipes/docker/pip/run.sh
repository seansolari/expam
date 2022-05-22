#!/bin/bash

# Modify source variable to where you would like
# database files to persist on your local computer.
docker run \
    --mount type=bind,source=/Users/ssol0002/Documents/Projects/pam/test/data,target=/app/INPUT \
    --mount type=bind,source="$(pwd)/OUTPUT",target=/app/OUTPUT \
    --shm-size 8G \
    -t expam-pip \
    \
