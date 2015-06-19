#!/usr/bin/env Rscript

#
# Reads mtx matrix txt file and row- and column names, respectively.
# Writes out a sparse matrix rds R object.
# 
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

# Load library
require(Matrix)

# Set paths
targetmatrix <- "/path/to/targetmatrix"
outpath <- "/path/to/rds

# Read matrix
X <- readMM(sprintf("%s.mtx", targetmatrix))

# Set row and col names
rn <- read.delim(sprintf("%s.rownames", targetmatrix), stringsAsFactors = FALSE, header = FALSE)
cn <- read.delim(sprintf("%s.colnames", targetmatrix), stringsAsFactors = FALSE, header = FALSE)
rownames(X)[rn[,1]] <- rn[,2]
colnames(X)[cn[,1]] <- cn[,2]

# Transform scores
X <-as(1 - 2^X, "dgTMatrix")
X[X < 0] <- 0

# Set on-target here. Information should be provided by your library vendor

# Remove 0 columns
X <- X[,which(colSums(X) != 0)]

# Write it out
saveRDS(X, file = outpath)
