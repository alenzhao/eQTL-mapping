######################
# Author: E. Schutte #
######################

### Settings
## Command line arguments
# To be implemented

## Loading packages
# DESeq used to process raw reads.
library("DESeq")

## Loading data
# Load data file with raw expression data.
load("~/Dropbox/Erik/expr_vst_condition_patient_rmBatch_88samples.Rdata")
class(expr.vst.rmBatch) <- "numeric"
cds = newCountDataSet(expr.vst.rmBatch, condition)
