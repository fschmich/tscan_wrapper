# Predict siRNA-to-gene off-target relation matrices using TargetScan

Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)

The following steps will guide you through a set scripts used for the prediction of siRNA off-targets and the construction of siRNA-to-gene target relation matrices as used for input to the gespeR model (www.cbg.ethz.ch/software/gespeR). It is required to install TargetScan and TargetScan context+ (http://www.targetscan.org/). Please adjust paths in each script according to your installation.

1. Prepare all seed sequences running the script "1_make_seeds.py".

   This script creates one file for each n-mer seed sequence, suitable for input to TargetScan's targetscan_60.pl input script

2. Predict all seed matches running script "2_run_targetscan_60.sh".

   This script loops over all seeds created in step 1. For consecutive predictions, i.e. for multiple libraries, step 1 and 2 only need to be done once. We have provided a file containing the 3' UTR sequences with with removed gaps for hg19 in the ./data folder.

3. Predict the context+ scores for all matches found in step 2 running script "3_run_predicttargets.py".

   This script should read your siRNA library file and loop over siRNA IDs and antisense sequences. As an example, we included three calls to siRNAs from Ambion.

4. Construct the siRNA-to-gene target relation matrix from all context+ score predictions in step 3 running script "4_run_makeX.sh".

   This script runs makeX.py and constructs a txt file for sparse matrices in the mtx format. As input, in addition to the context+ score predictions from step 3, it requires a file containing all siRNA IDs and all gene IDs, respectively. For the Ambion example, we provided the file with siRNA IDs and gene IDs in the ./data folder.

5. Construct a sparse matrix R object from the mtx txt file running script "5_mtx2rds.R". 

   This script takes as input the path/basename from the output of step 4 and reads the mtx file format. The script adds row- and column names and saves the final siRNA-to-gene target relations matrix as an rds file. We point out the section in the code, where the addition of the on-target component can be implemented.
