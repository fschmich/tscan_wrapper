#!/bin/sh

# Run TargetScan's base script (targetscan_60.pl) on all seed files.
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

SEEDS_DIR="/path/to/seeds"
TS_UTRSEQ="./data/hg19_ucsc_3p.txt"
TSCAN_60="/path/to/targetscan_60.pl"
OUT_DIR="/path/to/tscan60"

for s in $SEEDS_DIR/*.txt
do
	bname=$(basename $s)
	perl $TSCAN_60 $s $TS_UTRSEQ $OUT_DIR/tscan.$bname
done
