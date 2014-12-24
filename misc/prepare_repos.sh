#!/bin/bash
#to be start from root pkg folder 
cp DESCRIPTION PACKAGES
gzip PACKAGES
mv PACKAGES.gz ../www/R/src/contrib
cp inst/doc/GenometriCorr.pdf ../www
