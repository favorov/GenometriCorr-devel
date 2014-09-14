#!/bin/bash
cp DESCRIPTION PACKAGES
gzip PACKAGES
mv PACKAGES.gz ../www/R/src/contrib
