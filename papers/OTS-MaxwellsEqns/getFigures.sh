#!/bin/sh
#
# This script collects the figures from the src directory and places them
# in the figures directory.
#
# Kevin T. Chu
# February 2008
#

FIGURES_DIR="figures"

SRC_DIR="../../src"
EXAMPLE_DIR_LIST="1d/advection_diffusion_eqn                            \
                  1d/burgers_eqn/single_hump                            \
                  1d/diffusion_eqn                                      \
                  1d/fourth_order_parabolic                             \
                  2d/advection_eqn                                      \
                  2d/diffusion_eqn/square_domain                        \
                  2d/diffusion_eqn/irregular_domain/circle_domain       \
                  2d/diffusion_eqn/irregular_domain/starfish_domain     \
                 "

# create directory for figures if it doesn't already exist
if [ ! -d "${FIGURES_DIR}" ]; then
  mkdir -p ${FIGURES_DIR}
fi

# move figures 
for DIR in ${EXAMPLE_DIR_LIST}; do
  FIGURES=`ls ${SRC_DIR}/${DIR}/figures/*.eps 2> /dev/null` 
  if [ ! -z "${FIGURES}" ]; then
    echo -n "Moving figures from "
    echo "${SRC_DIR}/${DIR}/figures/ to ${FIGURES_DIR} ..."
    mv ${SRC_DIR}/${DIR}/figures/*.eps ${FIGURES_DIR}/
  else
    echo -n "Figures missing from "
    echo "${SRC_DIR}/${DIR}/figures/"
  fi
done

