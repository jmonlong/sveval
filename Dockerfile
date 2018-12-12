# Inspired by Bioconductor containers

FROM bioconductor/release_base2

MAINTAINER jmonlong@ucsc.edu

COPY . /tmp/sveval

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/sveval/

RUN R -f /tmp/sveval/install.R

RUN R -e "devtools::install('/tmp/sveval')"

CMD ["R"]