FROM rocker/r-base

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
                libcurl4-openssl-dev \
                libssl-dev 

RUN install2.r -e shiny shinydashboard data.table rhandsontable littler BiocManager RColorBrewer
RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r limma edgeR
RUN install2.r -e dplyr
RUN install2.r -e gplots
RUN install2.r -e ggplot2

RUN rm -rf /tmp/downloaded_packages

# Fix the time zone
RUN unlink /etc/localtime && \
    ln -s /usr/share/zoneinfo/Europe/Berlin /etc/localtime

RUN mkdir /root/RNAinspector
COPY app.R /root/RNAinspector
COPY vignette.html /root/RNAinspector
COPY tables/ /root/RNAinspector/tables/

EXPOSE 2525

CMD ["R", "-e", "shiny::runApp('/root/RNAinspector', host = '0.0.0.0', port = 2525)"]



