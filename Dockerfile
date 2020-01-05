FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    gcc \
    gfortran \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libxslt1-dev \
    python3-dev \
    zlib1g-dev

RUN curl https://bootstrap.pypa.io/get-pip.py -o - | python3

WORKDIR pythomics
COPY Makefile MANIFEST.in setup.py setup.cfg tox.ini ./
COPY requirements requirements
RUN pip3 install -r requirements/testing.txt -r requirements/linux.txt

COPY scripts scripts
COPY tests tests
COPY pythomics pythomics

RUN pip3 install -e .
