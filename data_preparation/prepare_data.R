######################
# Author: E. Schutte #
######################

### Settings
## Command line arguments
# Clear the current environment
rm(list=ls())

# Makes command line arguments availble for this script.
args = commandArgs(trailingOnly = T)

# Test if there is at least two arguments: if not, return an exception.
if ( length(args) <= 1 )  {
  stop("At least two argument must be supplied!\nRscript --vanilla prepare_data.R <genotype_file> <gene_expression_file>\n", call.=F)
} else {
  # Let the user know the program started.
  cat("Starting prepare_data.R..\n\nLoading files..\n\n")
  
  # Load user provided files.
  load(args[1])
  load(args[2])
  #load("~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata")
  #load("~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata")
  # Print if the files are loaded.
  cat("Files loaded:\n",args[1],"\n",args[2],"\n\n")
  
  # Store all the names in the global environment.
  nms <- ls(pattern = "", envir = .GlobalEnv)
  L <- sapply(nms, get, envir = .GlobalEnv, simplify = FALSE)
  
  # Sort out the matrixes.
  L <- L[which(sapply(L,is.matrix) == T)]
  GE <- L[[1]]
  snps <- L[[2]]
}
### Functions
## Change_allele_to_frequencies
# Function to identify the genotype and determine wheter it is major, minor or heterozygote.
change_allele_to_frequencies <- function (x,y) {
  if ( identical(x,allele.list[[1]]$genotype) ) {
    # The genotype is major, return 0.
    return(0)
    
  } else if ( identical(allele.list[[4]]$genotype,x) ) {
    # The genotype is NA or 00, return NA.
    return('NA')
    
  } else if ( identical(x,allele.list[[3]]$genotype) ) {
    # If the genotype is heterozygote, return 3.
    return(1)
    
  } else {
    # If the genotype is one of the minor alleles, return 2.
    return(2)
  }
}

### Running the data preperation
## Loading files
# Main directory, should be universal on every system.
# Verbose - Print for checking existing paths.
cat("\nCheck if directorys exist..\n\n")
mainDir <- "~/Documents"

# Sub directory.
subDir <- "R/"

# If the main and sub directory do not exist, create them.
if ( dir.exists(file.path(mainDir, subDir)) ) {
  cat("Directory \'", file.path(mainDir, subDir),"\' already exists.",sep = "")
} else {
  dir.create(file.path(mainDir, subDir))
  cat("Directory \'", file.path(mainDir, subDir),"\' created!",sep="")
}
## Transform the snps
# Verbose - Print transforming SNPs.
cat("Transforming SNPs\n\n")
# For readability change the matrix name and remove Norway data.
snps <- mergedGenotypes_clones2[-(23:27),]

# Keeps track of the snps that are causing noise in the data and are not representable.
snps.discarded <- c()
snps.discarded.pos <- c()
snps.discarded.counter <- 0

# For every column in the genotype data.
for (i in 1:ncol(snps)) {
  
  # Verbose - Print current column and column name.
  #cat("#",i,"- current column: ",colnames(snps)[i],"\n")
  
  # Save genotype and genotype count for each column.
  alleles.all <- table(snps[,i])
  
  # A minimum of 3 samples for a genotype is required to effectively calculate correlation.
  # Therefore we are removing the rows where this is not the case.
  if ( length(table(snps[,i])) < 3 & min(table(snps[,i])) < 3) {
  
    # Set an index for the SNP, so that later it can be removed.
    snps.discarded.pos <- c(snps.discarded.pos,i)
    
    # Store the discarded SNP for later review.
    snps.discarded <- c(snps.discarded,colnames(snps)[i])
    snps.discarded.counter <- snps.discarded.counter + 1
    
  }
  
  # List with all the major, minor, and heterozygote genotypes and their respecitve counts.
  allele.list <- list(major <- list(genotype ="",count=0),
                      minor <- list(genotype =c(),count=c()),
                      hetero <- list(genotype ="",count=0),
                      other <- list(genotype ="",count=0))
  
  # Temp value for checking which genotype occurs the most.
  highest_count_homozygote <- 0
  
  # For every allele that the genotype is counted for.
  for (a in 1:length(alleles.all)) {
    # Extract the genotype.
    geno <- names(alleles.all[a])
    
    # Verbose - Print current allele.
    #cat("Current allele = ", geno,'\n')
    
    # Extract the number of occurrences.
    count <- alleles.all[[a]]
    
    # Split the genotype and check if the first allele is the same as the second
    # thus checking for heterozygote or homozygote genotype.
    alleles <- strsplit(geno, "")
    
    if (identical(alleles[[1]][1],alleles[[1]][2]) & geno != '00') {
      # Verbose - Print if alleles are identical, and thus are homozygous.
      #cat('Alleles are identical, sorted as homozygote\n')
      
      if (highest_count_homozygote == 0) {
        # For the first iteration, the first genotype's occurences are the most occured
        # genotype and will be added on the first position.
        
        # Verbose - Print which genotype and count is added as major allele.
        #cat("Adding ",geno," with count ", count, "as major allele\n")
        highest_count_homozygote <- count
        allele.list[[1]]$genotype <- geno
        allele.list[[1]]$count <- count
        
      } else if (highest_count_homozygote < count) {
        # If the current genotype has more occurences than the lastly noted genotype,
        # it is automatically the major allele on the first position.
        
        # Verbose - Print which genotype and count is added as new major allele,
        # and which genotype and count previously occupied that spot.
        #cat("Adding ",geno," with count ", count, "as major allele\nPreviously ",
        #allele.list[[1]]$genotype," with count ", allele.list[[1]]$count,'\n')
        
        # Altering the position from the previous current major allele to minor allele
        # in the allele.list.
        allele.list[[2]]$genotype <- c(allele.list[[2]]$genotype, allele.list[[1]]$genotype)
        allele.list[[2]]$count <- c(allele.list[[2]]$count, allele.list[[1]]$count)
        
        # Setting the new major allele on the first position.
        highest_count_homozygote <- count
        allele.list[[1]]$genotype <- geno
        allele.list[[1]]$count <- count
        
      } else {
        # If the highest occuring genotype already found has more occurences,
        # the resulting genotype will be the minor allele which is the 2nd position
        # in the list.
        
        # Verbose - Print which genotype and count is added as minor allele.
        #cat("Adding ",geno," with count ", count, "as minor allele\n")
        allele.list[[2]]$genotype <- c(allele.list[[2]]$genotype, geno)
        allele.list[[2]]$count <- c(allele.list[[2]]$count, count)
      }
      
    } else {
      if (identical(geno,'00')) {
        # Separating the '00'/NA's from the data, these will be represented as 9.
        allele.list[[4]]$genotype <- geno
        allele.list[[4]]$count <- count
        
      } else {
        # If the genotype is not '00' it is automatically a heterozygote.
        
        # Verbose - Print if genes are not identical, and thus are heterozygote.
        #cat('alleles are not identical, heterozygote\n')
        
        # Verbose - Print which genotype and count is added as heterozygous.
        #cat("Adding ",geno," with count ", count, "as heterozygote allele\n")
        allele.list[[3]]$genotype <- geno
        allele.list[[3]]$count <- count 
      }
    }
  }
  
  # For every row in the genotype data.
  for (j in 1:nrow(snps)) {
    # Check the current snp agains the determined major, minor and heterozygote genotypes,
    # and determine the allele frequency.
    snp.freq <- change_allele_to_frequencies(snps[j,i],allele.list)
    
    # Verbose - Print the current snip' genotype and the corrosponding frequency.
    #cat("Changing current snip: ",snps[j,i]," to: ", snp.freq,"\n")
    
    # Save the SNP frequency at previous occupied genotype location.
    snps[j,i] <- snp.freq
  }
}
# Verbose - Print removing discarded SNPs
cat("Removing discarded SNPs\n\n")

# Remove discarded snps from data set.
snps <- snps[,-snps.discarded.pos]

# Convert genotype matrix to numeric values for Matrix eQTL analaysis.
class(snps) <- "numeric"

# Transfer the genotype matrix so that the columns 'align' with the gene expression matrix.
snps.t <- t(snps)

# Remove those samples that don't occur in the gene expression data. 
# “TCC-09-1" and ”TCC-08-1”.
indices.snps <- c()
for (colname in colnames(GE)){
  # Save indices, where the colnames from the gene expression co-exist in the 
  # genotype matrix.
  strs <- strsplit(colname,"__")
  indices.snps <- c(indices.snps,grep(strs[[1]][2],colnames(snps.t)))
}

# Re-arange data frame.
snps.t <- snps.t[,unique(indices.snps)]

## Order the gene expression data
# Save the indices for the swap.
indices.ge <- c()

# For every colname in colnames of the genotype samples.
for (colname in colnames(snps.t)) {
  # Save indices, where the colnames from the genotype matrix co-exist in the 
  # gene expression matrix.
  indices.ge <- c(indices.ge,  grep(colname, colnames(GE)))
}

# Re-arange data frame
GE <- GE[,indices.ge]

### Save image
## Save the current data imgae for further use.
# Save this particular data image for basic eqtl mapping.
#curDate <- date()
#save.image(file=paste(file.path(mainDir, subDir),curDate,".RData",sep=""))
#save(GE,snps.t, file=paste(file.path(mainDir, subDir),curDate,".RData",sep=""))


### User information
## Prints information to the command line for the user.
# Shows the user how many SNPs are discarded and which ones.
cat("\nThe following SNPs were discarded due to lack of significance in occurence;\n\nSNP Names;\n")
for (snp.disc in snps.discarded) {cat(snp.disc,"\n")}
cat("\n\nTotal number of SNPs discarded;\n",snps.discarded.counter,"\n\n")

### Settings
## Loading libraries
# Matrix eQTL mapping.
library("MatrixEQTL")

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
time_intervals <- list(t0 = seq(1,dim(GE)[2],4),
                       t1 = seq(2,dim(GE)[2],4),
                       t2 = seq(3,dim(GE)[2],4),
                       t3 = seq(4,dim(GE)[2],4))

# Change gene expression data to 4 segments for the time intervals.
for (interval in time_intervals) {
  # Pattern for output file name.
  pattern_name <- paste("Basic eQTL - Threshold ",pvOutputThreshold," ",date())
  
  # Output file name and location.
  output_file_name = tempfile(pattern = pattern_name,tmpdir="~/Documents/R");
  
  # Load gene expression data.
  gene = SlicedData$new();
  gene$CreateFromMatrix(GE[,interval])
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
  unlink(output_file_name)
  
  # Results.
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n\n\n')
  cat('Detected eQTLs:', '\n')
  show(me$all$eqtls)
  
  # Plot the histogram of all p-values.
  plot(me)
}

#plot(me$all$eqtls)
#plot(me$all$eqtls$pvalue,me$all$eqtls$snps)

#me$all$eqtls$beta_se = me$all$eqtls$beta / me$all$eqtls$statistic
#plot(me$all$eqtls$beta_se)

