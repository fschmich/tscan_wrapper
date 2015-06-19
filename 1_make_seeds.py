#!/usr/bin/env python

# Create all seeds of length SEED_LEN. These files will serve as input to TargetScan's base script (targetscan_60.pl).
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

# Libraries
from itertools import product

SPECIES = "9606"
SEED_LEN = 7
OUT_DIR = "/path/to/seeds"

alphabet = ["A", "C", "G", "U"]
seeds = [''.join(i) for i in product(alphabet, repeat=SEED_LEN)]
for s in seeds:
 f = open("%s/%s.txt" % (OUT_DIR, s), "w")
 f.write("%s\t%s\t%s\n" % (s, s, SPECIES))
 f.close()
