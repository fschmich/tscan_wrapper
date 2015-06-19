#!/bin/sh

# Run predicttargets.py on a three example Ambion siRNAs.
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

OUT_DIR="/path/to/targetpredict"

python predicttargets.py --id s223211 --seq AAAAUACAACUGGAGAUGGag --out $OUT_DIR
python predicttargets.py --id s6087 --seq UUUCAGAUCUCGGUAGACGgt --out $OUT_DIR
python predicttargets.py --id s3076 --seq UAUAGUCUUCUUUUAGUCCag --out $OUT_DIR