#!/bin/sh

cd "$(dirname "$0")"

# do data preparation for DESeq2 analysis in R
## If running on other pc - make sure to check curdir used by python script
## and make sure that all data files are in the dir as the bash script
python3 read_count_class.py -R

# do DESeq2 analysis using all volunteers
# non-responders are all assigned factors as if they have been measured on day 0
Rscript deseq2_NDN.R

rm Rplots.pdf
