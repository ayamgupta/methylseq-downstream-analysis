FROM rocker/r-ver:4.4.3 AS builder

# Install system dependencies needed for building R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libbz2-dev \
    libzstd-dev \
    liblzma-dev \
    libgit2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean

# Install R packages (CRAN + Bioconductor + GitHub)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_version('BiocManager', version = '1.30.26', repos = 'https://cloud.r-project.org')" \
 && R -e "remotes::install_version('devtools', version = '2.4.5', repos = 'https://cloud.r-project.org')" \
 && R -e "remotes::install_version('dplyr', version = '1.1.4', repos = 'https://cloud.r-project.org')" \
 && R -e "remotes::install_version('purrr', version = '1.0.4', repos = 'https://cloud.r-project.org')" \
 && R -e "remotes::install_version('tidyr', version = '1.3.1', repos = 'https://cloud.r-project.org')" \
 && R -e "remotes::install_version('readr', version = '2.1.5', repos = 'https://cloud.r-project.org')" \
 && R -e "remotes::install_version('argparse', version = '2.2.5', repos = 'https://cloud.r-project.org')" \
 && R -e "BiocManager::install(version='3.20', ask=FALSE, update=FALSE)" \
 && R -e "BiocManager::install(c('genomation', 'GenomicRanges', 'txdbmaker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'ChIPpeakAnno', 'org.Hs.eg.db', 'enrichR'), ask=FALSE, update=FALSE, clean=TRUE)" \
 && R -e "remotes::install_github('al2na/methylKit@ecf85842bdd5252161fa5bbe6aa0d285c87e1f77')"

# Stage 2: Clean runtime image
FROM rocker/r-ver:4.4.3

# Install only runtime libraries (not headers or dev tools)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libbz2-dev \
    libzstd-dev \
    liblzma-dev \
    libgit2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean

# Copy installed R packages from builder stage
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

WORKDIR /data



# Single Stage build Can be used if you want to keep it simple, but it will include build dependencies in the final image.

# FROM rocker/r-ver:4.4.3

# # Install system dependencies
# RUN apt-get update && apt-get install -y \
#     libcurl4-openssl-dev \
#     libssl-dev \
#     libxml2-dev \
#     libbz2-dev \
#     libzstd-dev \
#     liblzma-dev \
#     libgit2-dev \
#     libxt-dev \
#     libharfbuzz-dev \
#     libfribidi-dev \
#     libfreetype6-dev \
#     libpng-dev \
#     libtiff5-dev \
#     libjpeg-dev \
#     && apt-get clean

# # Install CRAN packages and BiocManager
# RUN R -e "install.packages(c('BiocManager', 'devtools', 'dplyr', 'purrr', 'tidyr', 'readr'  ), repos='https://cloud.r-project.org')"

# # Install Bioconductor packages
# RUN R -e "BiocManager::install(version='3.20', ask=FALSE, update=FALSE)" \
#     && R -e "BiocManager::install(c('remotes', 'genomation', 'GenomicRanges', 'txdbmaker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'ChIPpeakAnno', 'org.Hs.eg.db', 'enrichR'), ask=FALSE, update=FALSE, clean=TRUE)"

# # Install methylKit from GitHub (latest dev version)
# RUN R -e "remotes::install_github('al2na/methylKit')"

# WORKDIR /data