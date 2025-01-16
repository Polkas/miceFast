# devcontainer/Dockerfile
FROM rocker/r-ver:4.3.0  # or whichever version you need

# Install system libs needed by Rcpp packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo'))"
