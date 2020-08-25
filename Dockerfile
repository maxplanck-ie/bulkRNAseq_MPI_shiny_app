FROM centos:7.2.1511

RUN yum install -y epel-release && \
    yum install -y R-3.6.0-1.el7.x86_64 wget nano vim-X11 vim-common vim-enhanced vim-minimal ypbind yp-tools ypserv autofs nfs-utils rsyslog && \
    yum install -y openssl-devel curl libcurl-devel mesa-libGLU freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel cairo-devel libxml2-devel xorg-x11-server-devel libX11-devel libXt-devel  && \
    mkdir /data && \
    mkdir /etc/automount && \
    R -e "install.packages(c('shiny','crosstalk','rmarkdown'), repos='https://cran.rstudio.com/',dependencies=TRUE)" 

# Fix the time zone
RUN unlink /etc/localtime && \
    ln -s /usr/share/zoneinfo/Europe/Berlin /etc/localtime

ADD ./mounts.py /usr/local/bin/mounts.py
ADD ./startup.sh /usr/local/bin/startup.sh

VOLUME ["/export/"]

RUN mkdir /root/bulkRNAseq
#COPY ui.R /root/snakequest
#COPY server.R /root/snakequest
COPY app.R /root/bulkRNAseq

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 2525

CMD ["/usr/local/bin/startup.sh"]



