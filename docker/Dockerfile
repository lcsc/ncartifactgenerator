ARG VARIANT=4.3
ARG BASE_IMAGE=tidyverse
# ghcr.io/rocker-org/devcontainer/tidyverse:4.3
FROM ghcr.io/rocker-org/devcontainer/${BASE_IMAGE}:${VARIANT}

ARG HDF5_BRANCH=chunk_functions
ARG NCARTIFACTGENERATOR_BRANCH=master

COPY *.sh /usr/local/bin/Rinstall

RUN apt update && apt install -y patch libhdf5-dev

RUN pwd && mkdir /tmp/hdf5r && cd /tmp/hdf5r && \
    curl https://github.com/HDFGroup/hdf5/releases/download/hdf5-1_14_1-2/hdf5-1_14_1-2-ubuntu-2204.tar.gz -OL && \
    tar -xvzf hdf5-1_14_1-2-ubuntu-2204.tar.gz && \
    tar -xvzf hdf5/HDF5-1.14.2.1-Linux.tar.gz && \
    cp -av HDF5-1.14.2.1-Linux/HDF_Group/HDF5/1.14.2.1/* /usr/


RUN INSTALL_SYS_REQS="true" PAK_VERSION="auto" \
    PACKAGES="github::lcsc/hdf5r@${HDF5_BRANCH},sf,rworldmap,chron,bit64,unix,ggimage" \
    Rinstall

RUN INSTALL_SYS_REQS="true" PAK_VERSION="auto" \
    PACKAGES="github::lcsc/ncartifactgenerator@${NCARTIFACTGENERATOR_BRANCH}" \
    Rinstall