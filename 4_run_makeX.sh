#!/bin/sh

# Run makeX.py on the output directory of predicted context+ scores for siRNA off-targets.
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

CPS_DIR="/path/to/targetpredict"
SIRNAS="./data/sirnas.txt"
SIRNAS="./data/genes.txt"

OUT="/path/to/targetmatrix"

python makeX.py --cps $CPS_DIR --sirnas $SIRNAS --genes $GENES --out $OUT