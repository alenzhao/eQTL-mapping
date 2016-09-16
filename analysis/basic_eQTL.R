######################
# Author: E. Schutte #
######################

### Settings
## Loading libraries
# Matrix eQTL mapping.
library("MatrixEQTL")

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
pvOutputThreshold = 1e-5;

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
  pattern_name <- paste("Threshold ",pvOutputThreshold," ",date())
  
  # Output file name and location.
  output_file_name = tempfile(pattern = pattern_name,tmpdir="~/Documents/R-Output");
  
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
  me = Matrix_eQTL_engine(
    snps = snps.sd,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
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
