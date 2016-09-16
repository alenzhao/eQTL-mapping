######################
# Author: E. Schutte #
######################

### Settings
## Loading libraries
# Matrix eQTL mapping.
library("MatrixEQTL")

## Command line arguments
# Makes command line arguments availble for this script.
args = commandArgs(trailingOnly = T)

# Test if there is at least one argument: if not, return an exception.
if ( length(args) == 0 )  {
  stop("At least two argument must be supplied!\nRscript --vanilla <gene_positions_file> <snps_positions_file> <gene_expression_file>\n", call.=F)
} else if ( length(args) <= 2 ) {
  #default output
  args[3] = "out.txt"
} else {
  stop("This is for testing, I don't know what will happen")
}

## Load default settings
# Main directory, should be universal on every system.
mainDir <- "~/Documents/"

# Sub directory.
subDir <- "/R-Output/"

# If the main and sub directory do not exist, create them.
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# The basic_eqtl.RData file is stored in project root/data_preparation/RData/basic_eqtl_mapping.RData
load("data_preparation/RData/basic_eqtl_mapping.RData")

### Prepare matrix eqtl
## Settings
# Set the model used for eqtl mapping.
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Change covariates to character() for no covariates.
covariates_file_name = character()

# Only associations significant at this level will be saved.
pvOutputThreshold_tra = 1e-5;
pvOutputThreshold_cis = 1e-5;

# Distance for local gene-SNP pairs
cisDist = 1e6;

# Set gene and snp position files.
snpspos = read.table("~/Documents/CeD_43loci_alternative_format.txt", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("~/Dropbox/Erik/genepos.txt", header = TRUE, stringsAsFactors = FALSE);
#genepos <- args[1]
#genepos <- "~/Dropbox/Erik/genepos.txt"
#snpspos <- args[2]
#snpspos <- "~/Documents/CeD_43loci_alternative_format.txt"

# Error covariance matrix.
# Set to numeric() for identity.
errorCovariance = numeric();

# Load genotype data.
snps.sd = SlicedData$new()
snps.sd$CreateFromMatrix(snps.t)
snps.sd$fileDelimiter = "\t";      # the TAB character
snps.sd$fileOmitCharacters = "NA"; # denote missing values;
snps.sd$fileSkipRows = 1;          # one row of column labels
snps.sd$fileSkipColumns = 1;       # one column of row labels
snps.sd$fileSliceSize = 2000;      # read file in slices of 2,000 rows

# Analysis genotype data versus gene expression data requries a loop through the 
# different timepoints measured in the gene expression data.
time_intervals <- list(t0 = seq(1,dim(count.all)[2], 4),
                       t1 = seq(2,dim(count.all)[2],4),
                       t2 = seq(3,dim(count.all)[2],4),
                       t3 = seq(4,dim(count.all)[2],4))

# Change gene expression data to 4 segments for the time intervals.
for (interval in time_intervals) {
  # Pattern for output file name.
  pattern_name.cis <- paste("Cis eQTL - Threshold ",pvOutputThreshold_cis," ",date())
  pattern_name.trans <- paste("Trans eQTL - Threshold ",pvOutputThreshold_tra," ",date())
  
  # Output file name and location.
  output_file_name_cis = tempfile(pattern = pattern_name.cis,tmpdir="~/Documents/R-Output/cis");
  output_file_name_tra = tempfile(pattern = pattern_name.trans,tmpdir="~/Documents/R-Output/trans");
  
  # Load gene expression data .
  gene = SlicedData$new();
  gene$CreateFromMatrix(count.all[,interval])
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  
  # Load covariates.
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  
  ## Run the analysis.
  me = Matrix_eQTL_main(
    snps = snps.sd,
    gene = gene, 
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  # Unlink output_file_name (keep it commented to keep eQTL mappings),
  #unlink(output_file_name)
  
  # Results.
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n\n\n')
  #cat('Detected eQTLs:', '\n')
  #show(me$all$eqtls)
  
  # Plot the histogram of all p-values.
  plot(me)
}

plot(me$all$eqtls)
plot(me$all$eqtls$pvalue,me$all$eqtls$snps)
