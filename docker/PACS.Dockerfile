FROM bioconductor/bioconductor_docker:3.18-R-4.3.3
RUN Rscript -e "install.packages( \
    c('devtools'),\
    repos = BiocManager::repositories())"
RUN Rscript -e "devtools::install_github('Zhen-Miao/PICsnATAC', repos = BiocManager::repositories())"
RUN Rscript -e "devtools::install_github('Zhen-Miao/PACS', repos = BiocManager::repositories())"
