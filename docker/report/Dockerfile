FROM rocker/r-base:4.4.2

COPY install.R /scripts/

RUN apt update && \
    mkdir -p /install && \
    wget -O /install/quarto-1.8.27-linux-amd64.deb https://github.com/quarto-dev/quarto-cli/releases/download/v1.8.27/quarto-1.8.27-linux-amd64.deb && \
    dpkg -i /install/quarto-1.8.27-linux-amd64.deb && \
    apt install -y procps \
        python3-yaml \
        libcurl4-openssl-dev \
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
