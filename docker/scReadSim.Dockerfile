FROM python:3.10-slim
RUN apt-get update && apt-get install -y procps
RUN apt-get install gcc
RUN pip install --no-cache-dir cykhash
RUN apt install zlib1g-dev
RUN pip install macs3
RUN apt install r-base
RUN apt install wget
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 -O samtools.tar.bz2
RUN tar -xjvf samtools.tar.bz2
RUN apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
RUN tar -zxvf bedtools-2.29.1.tar.gz
RUN cd bedtools2
RUN make
RUN apt install libxml2-dev
RUN apt install fontconfig libfontconfig1-dev libfreetype6 libfreetype6-dev
RUN apt install libharfbuzz-dev libfribidi-dev
RUN apt install libtiff5-dev
