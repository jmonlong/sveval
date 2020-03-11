FROM ubuntu:18.04

MAINTAINER jmonlong@ucsc.edu

## Dependencies
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        wget \
        screen \
        less \
        nano \
        git \
        curl \
        python3 \
        python3-dev \
        python3-pip \
        python3-setuptools \
        gcc \
        make \
        bzip2 \
        tabix \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libxml2-dev \
        libssl-dev \
        libmariadbclient-dev \
        libcurl4-openssl-dev \
        apt-transport-https \
        software-properties-common \
        dirmngr \
        gpg-agent \
        && rm -rf /var/lib/apt/lists/*

## Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
        && apt-get update \
        && DEBIAN_FRONTEND=noninteractive apt-get install -y r-base r-base-dev

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/sveval/

## Install sveval
COPY . /tmp/sveval
RUN R -e "install.packages('BiocManager')" \
        && R -f /tmp/sveval/install.R \
        && R -e "devtools::install('/tmp/sveval')"

## Install bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
        tar -xjf bcftools-1.10.2.tar.bz2 && \
        cd bcftools-1.10.2 && \
        ./configure && make && make install && \
        cd .. && rm -rf bcftools-1.10.2 bcftools-1.10.2.tar.bz2

## Install snakemake
RUN pip3 install --upgrade pip
RUN pip3 install --no-cache-dir snakemake==5.8.2

WORKDIR /home
