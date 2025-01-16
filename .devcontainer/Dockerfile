FROM rocker/r-ver:4.3.0

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    libgit2-dev \
    build-essential \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'devtools', 'data.table', 'dplyr', 'ggplot2', 'mice', 'rmarkdown', 'knitr'))"
