#!/bin/sh
#
# This script generates PDF figures from the eps versions of figures.
#
# Kevin T. Chu
#

FIGURES_DIR="figures"

# convert figures
cd ${FIGURES_DIR}
for FIG in `ls *.eps`; do
  echo -n "Converting ${FIG} ..."
  epstopdf ${FIG}
  echo "done"
done

