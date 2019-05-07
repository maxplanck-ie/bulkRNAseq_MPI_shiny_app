#!/bin/bash

#Set up the mounting
python /usr/local/bin/mounts.py

#R --no-save --no-restore
#chown -R galaxy:bioinfo /root/bulkRNAseq
#su galaxy -s /bin/bash -c 'R -e ''shiny::runApp\(\"/root/bulkRNAseq\",port=2525,host=\"0.0.0.0\"\)''' 
R -e 'shiny::runApp("/root/bulkRNAseq",port=2525,host="0.0.0.0")'
