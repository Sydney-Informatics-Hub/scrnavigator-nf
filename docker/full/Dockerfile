FROM rocker/r-base:4.4.2

COPY install.R /scripts/

RUN apt update && \
    apt install -y libcurl4-openssl-dev \
        libssl-dev \
        libfontconfig1-dev \
        libxml2-dev \
        libhdf5-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libwebp-dev \
        libv8-dev && \
    Rscript --vanilla /scripts/install.R