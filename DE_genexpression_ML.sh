#!/bin/sh

cd "$(dirname "$0")"

echo "Running ML based pipeline"

# do data preparation for DESeq2 analysis in R
python3 readcountclass.py -ML

# do DESeq2 analysis using all volunteers
# non-responders are all assigned factors as if they have been measured on day 0
#Rscript deseq2_NDN.R

