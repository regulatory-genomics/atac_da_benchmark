# Build image
# Use an official Ubuntu base image
FROM ubuntu:20.04

# Set non-interactive installation to avoid getting stuck with prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary tools and libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gpg-agent \
    lsb-release \
    build-essential \
    wget \
    curl \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    python3 \
    python3-pip \
    python3-dev

# Add the CRAN repository for the latest R versions
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# Install R
RUN apt-get update && \
    apt-get install -y --no-install-recommends r-base

# Clean up to reduce image size
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory (optional)
WORKDIR /workspace

# Default command: open a bash shell (change as needed)
CMD ["/bin/bash"]

libfontconfig1-dev
apt-get install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev 
libharfbuzz-dev libfribidi-dev
apt-get install gfortran
apt-get install libblas-dev liblapack-dev
RUN apt-get update && apt-get install -y \
    procps \
devtools::install_github("SONGDONGYUAN1994/scDesign3")
