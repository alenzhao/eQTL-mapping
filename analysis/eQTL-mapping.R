######################
# Author: E. Schutte #
######################

# Let the user know the packages are loaded
cat("Loading packages, if not present, they are installed.\n\nThis may take a while.\n\n")
### Libraries
## Source.
# Set source(s) needed.
source("http://bioconductor.org/biocLite.R")

## CRAN
# Check if all libraries are installed, if not install them.
packages <- c("MatrixEQTL","devtools")
if ( length( setdiff( packages, rownames( installed.packages() ) ) ) > 0) {
  suppressMessages(
  install.packages( setdiff( packages, rownames( installed.packages() ) ) )  
  )
}

## BiocLite
# Load BiocLite packages.
bioclite.packages <- c("vegan", "Rsamtools", "qvalue")
if( length( setdiff( bioclite.packages, rownames( installed.packages() ) ) > 0 ) ) {
  suppressMessages(
  biocLite( setdiff( bioclite.packages, rownames( installed.packages() ) ) )
  )
}

## Libraries
# Load Matrix eQTL library.
library("MatrixEQTL")
library("devtools")
if ( length( dev_packages() == 0 ) ) {
  devtools::install_github("guigolab/sQTLseekeR")
} 
library("sQTLseekeR")

### Functions
## Global settings.
# Stores settings globally.
global <- function(x,y) {
  assign(x,y,envir=.GlobalEnv)
}
## Settings.
# Set all needed settings.
settings <- function() {
  
  ## Command line arguments.
  # Clear the current environment.
  rm(list=ls())
  
  # Makes command line arguments availble for this script.
  global("args",commandArgs(trailingOnly = T))
  
  # Command line Test 1:
  #global("args",list("-h"))
  
  # Command line Test 2:
  #global("args",list("basic","~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata"))
  
  # Command line Test 3:
  #global("args", list("basic","~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata","-v"))
  
  # Command line Test 4:
  #global("args", list("basic","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata","-v"))
  
  # Command line Test 5:
  #global("args",list("basic","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata"))
  
  # Command line Test 6:
  #global("args",list("~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata","-v"))
  
  # Command line Test 7:
  #global("args", list("~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata"))
  
  # Command line Test 8:
  #global("args",list("basic","~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata","~/Dropbox/Erik/expr_vst_condition_patient_rmBatch_88samples.Rdata"))
  
  #load("~/Dropbox/Erik/clones_CeD_genotypes_noDuplicates.Rdata")
  #load("~/Dropbox/Erik/count_all_batch_times_sample.names_96Samples_noCD8_leiden.Rdata")
  #load("~/Dropbox/Erik/expr_vst_condition_patient_rmBatch_88samples.Rdata")
  
  # Set mapping flag.
  global("MAPPING_TYPE", "")
  global("MAPPING_OPTS", c("basic","ct","s"))
  
  # Set verbose flag.
  global("VERBOSE", FALSE)
}

## Error handling.
# Handles error messages, exits script.
error <- function(error,input) {
  tmp <- ""
  if ( length(input) >= 2 ) {
    for ( i in 1:length(input) ) {
      tmp <- paste(as.character(input[[i]]), sep="")
    }
  }
  input <- tmp
  cat("
      Something went wrong:
      
      Your input:",
      input,"is probably invalid.
      Reason:",
      error,
      
      "
      Please check for spelling errors, and if the correct file/option is entered as input.
      ")
  q(save="no")
}

## Help.
# display help message.
help <- function() {
  cat("
  QTL-mapping Script
  
  This script maps genotypes from samples to expression values, locating basic-, cis-, trans-
  and s-QTLs. Important to note; The dimensions from the genotype file and the gene expression
  file have to be from the same length, or a multiplication of that same length (given time
  points).
  
  !!The data has to be saved in an *.Rdata image, that is then loaded into the script.!!
  
  Arguments:
  --help, -h                   - Displays this help message.
  
  map_type                     - The sort of eQTL Analysis you want performed
  |- Options:
  |- 'basic', basic eQTL analysis.
  |- 'ct', cis- trans-eQTL analysis.
  |- 's', sQTL analysis (NOT IMPLMENTED).
  
  genotype_file                - A File containing the genotypes for merged alleles per sample (column),
                                 per gene (row).
  
  gene_expression_file         - A File containing the gene expressino data per sample (column)
                                 per gene (row).
  
  --verbose, -v                - A Verbose option. Shows output of steps taken. You might want to save
                                 the std.out to a file ( see example ).
  
  Examples:
  - Texual;
  |-Rscript --vanilla eQTL-mapping.R <map_type> <genotype_file> <gene_expression_file> [--verbose, -v]
  
  - Script call;
  |-Rscript --vanilla eQTL-mapping.R basic myGenotypeFile.Rdata myGeneExpressionFile.Rdata
  |
  |- With verbose output;
  |  |-Rscript --vanilla eQTL-mapping.R basic myGenotypeFile.Rdata myGeneExpressionFile.Rdata -v
  |
  |- With verbose output, saved to a file;
  |-Rscript --vanilla eQTL-mapping.R basic myGenotypeFile.Rdata myGeneExpressionFile.Rdata -v > output.txt
  
  ")
  q(save="no")
}

## parseOpts.
# Parse command line options and test them.
parseOpts <- function() {
  
  # If length of arguments is 1.
  if ( length(args) == 1 ) {
    
    # Check if the first argument is a help flag.
    if ( length(args[1][[1]] != 0 | !is.null(args[1][[1]]) && args[1][[1]] == '-h' | args[1][[1]] == '--help' ) ) {
      help()
    # If it is not a help flag, throw error unknown flag.
    } else {
      error( "Unknown flag", args[1][[1]] )
    }
  
  # If length of arguments is 2 or 3.
  } else if ( length(args) == 2 | length (args) == 3 ) {
    
    # If the length of the first argument is not 0 or the first argument is not null.
    if ( length(args[1][[1]]) != 0 | !is.null(args[1][[1]]) ) {
      
      # If the first argument is in MAPPING_OPTS.
      if ( args[1][[1]] %in% MAPPING_OPTS ) {
        
        # Set first argument as global variable MAPPTING_TYPE.
        global( "MAPPING_TYPE", as.character(args[1][[1]]) )
        
      # If the first argument is not in MAPPING_OPTS throw error invalid mapping type.
      } else {
        error( "Invalid mapping type.",args[1][[1]] )
      }
    }
    
    # If argument 2 is zero or null and argument 3 is zero or null.
    if ( length(args[2][[1]]) != 0  && length(args[3][[1]]) != 0  ) {
      #| !is.null(args[2][[1]])  | !is.null(args[3][[1]])
      # Check if file argument 2 and argument 3 exist.
      if ( file.exists( file.path(args[2][[1]]) ) && file.exists( file.path(args[3][[1]]) ) ) {
        
        # Create temp space.
        temp.space <- new.env()
        
        # Load data into temp space.
        snps.temp <- load(args[2][[1]], temp.space)
        GE.temp <- load(args[3][[1]], temp.space)
        
        # Save data to current environment.
        snps <- get(snps.temp, temp.space) # Genotype file.
        GE <- get(GE.temp, temp.space) # Gene expression file.
        
        # Save data globally.
        global("GE", GE)
        global("snps", snps)
        
        # Store all the names in the global environment.
        #nms <- ls(pattern = "", envir = temp.space)
        #L <- sapply(nms, get, envir = .GlobalEnv, simplify = FALSE)
        
        # Sort out the matrixes.
        #L <- as.matrix(L[which(sapply(L,is.matrix) == T)])
        
      # If files don't exist, throw error.
      } else {
        
        # If the -v, --verbose flag is used as file input, throw an error.
        if ( args[2][[1]] == '-v' | args[2][[1]] == '--verbose' | args[3][[1]] == '-v' | args[3][[1]] == '--verbose' ) {
          error("-v, --verbose flag used as file input",args[2:3])
        }
        
        error("File/files does/do not exist",args[2:3])
      }
      
    # If argument 2 and 3 are empty, throw error to provide files.
    } else {
      error("Provide two input files",args[2:3])
    }
  
  # If length of args is 4.  
  } else if ( length (args) == 4 ) {
    
    # If the length of the first argument is not null or the first argument is not null.
    if ( length(args[1][[1]]) != 0 | !is.null(args[1][[1]]) ) {
      
      # If the first argument is in MAPPING_OPTS.
      if (args[1][[1]] %in% MAPPING_OPTS) {
        
        # Set the first argument as the MAPPING_TYPE.
        global("MAPPING_TYPE", as.character(args[1][[1]]))
      
      # If the first argument is not in MAPPING_OPTS, throw error invalid mapping type.  
      } else {
        error("Invalid mapping type.",args[1][[1]])
      }
    }
    
    # If the length of the 2nd argument is not null or the 2nd argument is not null,
    # and if the length of the 3th argument is not null or the 3th argument is not null.
    if ( length(args[2][[1]]) != 0  && length(args[3][[1]]) != 0  ) {
      
      # If the file as 2nd argument exists and the file as 3th argument exists.
      if( file.exists( file.path(args[2][[1]]) ) && file.exists( file.path(args[3][[1]]) ) ) {
        
        # Load data into environment
        temp.space <- new.env()
        GE <- get(load(args[2][[1]], temp.space)) # Genotype file.
        snps <- get(which(sapply(load(args[3][[1]], temp.space),is.matrix))) # Gene expression file.
        print(snps)
        
        
        # Store all the names in the global environment.
        #nms <- ls(pattern = "", envir = temp.space)
        #L <- sapply(nms, get, envir = .GlobalEnv, simplify = FALSE)
        
        # Sort out the matrixes.
        #L <- as.matrix(L[which(sapply(L,is.matrix) == T)])
        #global("GE",L[1])
        #global("snps", L[2])
        
      # If files don't exist, throw error.
      } else {
        
        # If the -v, --verbose flag is used as file input, throw an error.
        if ( args[2][[1]] == '-v' | args[2][[1]] == '--verbose' | args[3][[1]] == '-v' | args[3][[1]] == '--verbose' ) {
          error("-v, --verbose flag used as file input",args[2:3])
        }
        error("File/files does/do not exist",args[2:3])
      }
    
    # If the length of the 2nd argument is null or the 2nd argument is null,
    # and if the length of the 3th argument is null or the 3th argument is null,
    # throw an error to provide two input files.
    } else {
      error("Provide two input files",args[2:3])
    }
    
    # If the length of the 4th argument is not null or the 4th argument is not null.
    if ( length(args[4][[1]]) != 0 | !is.null(args[4][[1]]) ) {
      
      # If the 4th argument equals -v or --verbose.
      if ( args[4][[1]] == "-v" | args[4][[1]] == "--verbose" ) {
        
        # Set the VERBOSE variable to TRUE.
        global("VERBOSE", TRUE)
      
      # Throw an error of an unkown flag.  
      } else {
        error("Unknown flag",args[4][[1]])
      }
    }
  
  # If no arguments are provided throw an error.
  } else {
    error("No arguments provided",args)
  }
}

## Change_allele_to_frequencies.
# Function to identify the genotype and determine wheter it is major, minor or heterozygote.
change_allele_to_frequencies <- function (x,allele.list) {
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

## Data preparation.
data_prep <- function () {
  ### Data-preparation
  ## Transform the snps.
  # Verbose.
  if ( VERBOSE == TRUE ) {
    # Verbose - Print transforming SNPs.
    cat("Transforming SNPs\n\n")
  }
  
  # For readability change the matrix name and remove Norway data.
  snps <- snps[-(23:27),]
  
  # Keeps track of the snps that are causing noise in the data and are not representable.
  snps.discarded <- c()
  snps.discarded.pos <- c()
  snps.discarded.counter <- 0
  
  # For every column in the genotype data.
  for (i in 1:ncol(snps)) {
    # Verbose
    if ( VERBOSE == TRUE ) {
      # Verbose - Print current column and column name.
      cat("#",i,"- current column: ",colnames(snps)[i],"\n")
    }
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
      
      # Verbose
      if ( VERBOSE == TRUE ) {
        # Verbose - Print current allele.
        cat("Current allele = ", geno,'\n')
      }
      
      
      # Extract the number of occurrences.
      count <- alleles.all[[a]]
      
      # Split the genotype and check if the first allele is the same as the second
      # thus checking for heterozygote or homozygote genotype.
      alleles <- strsplit(geno, "")
      
      if (identical(alleles[[1]][1],alleles[[1]][2]) & geno != '00') {
        # Verbose
        if ( VERBOSE == TRUE ) {
          # Verbose - Print if alleles are identical, and thus are homozygous.
          cat('Alleles are identical, sorted as homozygote\n')
        }
        
        if (highest_count_homozygote == 0) {
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as major allele.
            cat("Adding ",geno," with count ", count, "as major allele\n")
          }
          
          # For the first iteration, the first genotype's occurences are the most occured
          # genotype and will be added on the first position.
          highest_count_homozygote <- count
          allele.list[[1]]$genotype <- geno
          allele.list[[1]]$count <- count
          
        } else if (highest_count_homozygote < count) {
          # If the current genotype has more occurences than the lastly noted genotype,
          # it is automatically the major allele on the first position.
          
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as new major allele,
            # and which genotype and count previously occupied that spot.
            cat("Adding ",geno," with count ", count, "as major allele\nPreviously ",
                allele.list[[1]]$genotype," with count ", allele.list[[1]]$count,'\n')
          }
          
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
          allele.list[[2]]$genotype <- c(allele.list[[2]]$genotype, geno)
          allele.list[[2]]$count <- c(allele.list[[2]]$count, count)
          
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as minor allele.
            cat("Adding ",geno," with count ", count, "as minor allele\n")
          }
        }
        
      } else {
        if (identical(geno,'00')) {
          # Separating the '00'/NA's from the data, these will be represented as 9.
          allele.list[[4]]$genotype <- geno
          allele.list[[4]]$count <- count
          
        } else {
          # If the genotype is not '00' it is automatically a heterozygote.
          allele.list[[3]]$genotype <- geno
          allele.list[[3]]$count <- count 
          
          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print if genes are not identical, and thus are heterozygote.
            cat('alleles are not identical, heterozygote\n')
          }
          
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as heterozygous.
            cat("Adding ",geno," with count ", count, "as heterozygote allele\n")
          }
        }
      }
    }
    
    # For every row in the genotype data.
    for (j in 1:nrow(snps)) {
      # Check the current snp agains the determined major, minor and heterozygote genotypes,
      # and determine the allele frequency.
      snp.freq <- change_allele_to_frequencies(snps[j,i],allele.list)
      
      # Verbose.
      if ( VERBOSE == TRUE ) {
        # Verbose - Print the current snip' genotype and the corrosponding frequency.
        cat("Changing current snip: ",snps[j,i]," to: ", snp.freq,"\n")
      }
      
      # Save the SNP frequency at previous occupied genotype location.
      snps[j,i] <- snp.freq
    }
  }
  
  # Verbose.
  if ( VERBOSE == TRUE ) {
    # Verbose - Print removing discarded SNPs
    cat("\n\nRemoving discarded SNPs:\n")
  }
  
  # Remove discarded snps from data set.
  snps <- snps[,-snps.discarded.pos]
  
  # Convert genotype matrix to numeric values for Matrix eQTL analaysis.
  class(snps) <- "numeric"
  
  # Transfer the genotype matrix so that the columns 'align' with the gene expression matrix.
  snps.t <- t(snps)
  
  # Change the sample names to the corrosponding format for the Gene Expression file.
  if ( length(grep("\\.",colnames((GE))) != 0) ) {
    colnames(GE) <- gsub("\\.", "_", colnames(GE))
  } else if ( length(grep("\\-",colnames((GE))) != 0) ) {
    colnames(GE) <- gsub("\\-", "_", colnames(GE))
  }
  
  # Change the sample names to the corrosponding format for the Genotype file.
  if ( length(grep("-",colnames((snps.t))) != 0) ) {
    colnames(snps.t) <- gsub("\\-", "_", colnames(snps.t))
  }
  
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
  
  # Save data.
  global("snps.t",snps.t)
  
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
  
  # Save data.
  global("GE",GE)
  
  ### User information
  ## Prints information to the command line for the user.
  # Shows the user how many SNPs are discarded and which ones.
  cat("\nThe following SNPs were discarded due to lack of significance in occurence;\n\nSNP Names;\n")
  for (snp.disc in snps.discarded) {cat(snp.disc,"\n")}
  cat("\nTotal number of SNPs discarded;\n",snps.discarded.counter,"\n\n")
}

## eQTL.
# Maps basic eQTLs.
eQTL <- function(){
  
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
                         t10 = seq(2,dim(GE)[2],4),
                         t30 = seq(3,dim(GE)[2],4),
                         t180 = seq(4,dim(GE)[2],4))
  
  # Empty eqtls variable.
  eqtls <- NULL
  
  # Change gene expression data to 4 segments for the time intervals.
  for (interval in 1:length(time_intervals) ) {
    
    # Pattern for output file name.
    pattern_name <- paste(mainDir,subDir,"lncrna_eqtl_basic_",names(time_intervals[interval])[1],".csv",sep="")
    file.create(pattern_name)
    
    # Output file name and location.
    output_file_name = pattern_name
    
    # Load gene expression data.
    gene = SlicedData$new();
    gene$CreateFromMatrix(GE[,time_intervals[[interval]]])
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
    
    # Verbose - Show all detected eQTLs.
    if ( VERBOSE == TRUE ) {
      cat('Detected eQTLs:', '\n')
      show(me$all$eqtls)
    }
    
    eqtls <- rbind(eqtls, me$all$eqtls)
    # Plot the histogram of all p-values.
    #plot(me)
  }
  
  #eqtls.rearanged <- cbind(eqtls[,4],eqtls[,1],eqtls[,2],eqtls[,5])
  #colnames(eqtls.rearanged) <- c("pvalue","SnpRs","EnsembleGeneID","FDR")
  #rownames(eqtls.rearanged) <- rownames(eqtls)
  eqtls <- format(eqtls,scientific=TRUE)
  
  cat("Writing to file...\n\n")
  write.csv(eqtls,file=paste(file.path(mainDir,subDir,"Basic",fsep=""), "/lncrna_eqtls_basic",sep=""))
  
  cat("Wrote file to: ",paste(file.path(mainDir,subDir,"Basic",fsep=""), "/lncrna_eqtls_basic",sep=""),"\n\n")
  cat("FINISHED\n\n")
  ### Login
  ## Login to Molgenis.
  # Login to retrieve the gene names.
  #source("http://localhost:8080/molgenis.R")
  #molgenis.login("admin", "admin")
  
  ## Prepared Statement.
  # Create query for data filtration.
  #qs <- ""
  #geneInformation = NULL
  #for (i in 1:length(eqtls.rearanged[,3])){
  #  qs[i] <- paste("EnsemblGeneID==", eqtls.rearanged[,3][i], sep="")
  #  if (i %% 500 == 0) {
  #    geneInformation <- rbind(geneInformation, molgenis.get("lncrna_GeneInfo", num=100000, q=q))
  #  }
  
  #}
  #q=paste(qs, collapse=",")
  
  ## Get data
  # Load the data from the Molgenis Repository using the Prepared statement.
  #geneInformation <- molgenis.get("lncrna_GeneInfo", num=100000, q=q)
}

## ctQTL.
# Maps cis- and trans-eQTLs
ctQTL <- function () {
  
  # Set gene and snp position files.
  snpspos = read.table("~/tijdelijk/eQTL-mapping/analysis/snp.txt", header = TRUE, stringsAsFactors = FALSE)
  genepos = read.table("~/tijdelijk/eQTL-mapping/analysis/gene.txt", header = TRUE, stringsAsFactors = FALSE)
  
  # Only associations significant at this level will be saved.
  pvOutputThreshold_tra = 1e-5
  pvOutputThreshold_cis = 0.05
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6
  
  ### Prepare matrix eqtl
  ## Settings
  # Set the model used for eqtl mapping.
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Change covariates to character() for no covariates.
  covariates_file_name = character()
  
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
                         t10 = seq(2,dim(GE)[2],4),
                         t30 = seq(3,dim(GE)[2],4),
                         t180 = seq(4,dim(GE)[2],4))
  
  # Empty variable for cis- trans- qtls.
  transqtls <- NULL
  cisqtls <- NULL
  
  # Change gene expression data to 4 segments for the time intervals.
  for (interval in time_intervals) {
    
    # Pattern for output file name.
    pattern_name.cis <- paste(mainDir,subDir,"lncrna_eqtl_cis_",names(time_intervals[interval])[1],".csv",sep="")
    pattern_name.trans <- paste(mainDir,subDir,"lncrna_eqtl_trans_",names(time_intervals[interval])[1],".csv",sep="")
    
    # Create output files.
    file.create(pattern_name.cis)
    file.create(pattern_name.trans)
    
    # Output file name and location.
    output_file_name_cis = pattern_name.cis
    output_file_name_tra = pattern_name.trans
    
    # Load gene expression data .
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
    unlink(output_file_name_cis)
    unlink(output_file_name_tra)
    
    # Results.
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n\n\n')
    
    # Verbose - Show all detected eQTLs.
    if ( VERBOSE == TRUE ) {
      cat('Detected eQTLs:', '\n')
      show(me$all$eqtls)
    }
    
    # Bind all cis- trans- eqtls to dataframe
    transqtls <- rbind(transqtls, me$trans$eqtls)
    cisqtls <- rbind(cisqtls, me$cis$eqtls)
    # Plot the histogram of all p-values.
    #plot(me)
  }
  transqtls <- format(transqtls,scientific=TRUE)
  cisqtls <- format(cisqtls,scientific=TRUE)
  
  cat("Writing to file...\n\n")
  write.csv( transqtls,file=paste(file.path(mainDir, subDir, "Trans", fsep = ""), "/lncrna_eqtls_trans", sep="") )
  write.csv( cisqtls,file=paste(file.path(mainDir, subDir, "Cis", fsep = ""), "/lncrna_eqtls_cis", sep="") )
  
  cat("Wrote file(s) to: \n",paste(file.path(mainDir, subDir, "Trans", fsep=""), sep=""),
      paste(file.path(mainDir, subDir, "Cis", fsep="") ) )
  
  cat("FINISHED\n\n")
}

## sQTL.
# Maps splice eQTLs.
sQTL <- function () {
  ### Pre-processing
  ## Change working directory.
  # Set working directory to test data location.
  setwd("/usr/local/lib/R/3.3/site-library/sQTLseekeR/data/")
  
  ## Load data.
  # Transcript expression for 5 genes and SNPs in these 5 regions, the third file contains the coordinates.
  trans.exp.f = "transExpression.tsv.gz"
  gene.bed.f = "genes.bed"
  genotype.f = "snps-012coded.tsv"
  
  ## Order data.
  # Ordered genotype file should be compressed and indexed if not done before.
  genotype.indexed.f = index.genotype(genotype.f)
  
  ## Read data.
  # Import transcript expression, cleaned.
  te.df = read.table(trans.exp.f, as.is=T, header=T, sep="\t")
  tre.df = prepare.trans.exp(te.df)
  
  # Show data.
  tre.df[1:5,1:5]
  
  ### Test for gene/SNP associations
  ## Tests.
  # Run test with transcript expression and genotype file and gene coordinates.
  gene.bed = read.table(gene.bed.f, as.is=T, sep="\t")
  
  # Set column names.
  colnames(gene.bed) = c("chr","start","end","geneId")
  
  res.df = sqtl.seeker(tre.df, genotype.indexed.f, gene.bed, svQTL=T, verbose=F)
  head(res.df)
  
  write.table(res.df, file="~/Documents/R/sQTLs-all.tsv", quote=FALSE, row.names=FALSE, sep="\t")
  
  sqtls.df = sqtls(res.df = res.df, FDR = 0.01, out.pdf="~/Documents/R/sQTLs-FDR01.pdf")
  head(sqtls.df)
}
  
## createDir.
# creates direcotry for mapping types.
createDir <- function(dir) {
  
  # If directory exists..
  if ( dir.exists( file.path(mainDir, subDir, dir, fsep="") ) ) {
    
    # .. show the user the directory already exists.
    cat("Directory: \'", file.path(mainDir, subDir, dir, fsep=""), " exists.\n", sep = "") 
  } 
  
  # If direcotry does not exists..
  else {
    
    # .. show the user the directory is being created.
    cat("Creating directory: \'", file.path(mainDir, subDir, dir, fsep=""),"\'\n",sep="")
    
    # Create the directory.
    dir.create( file.path(mainDir, subDir, dir, fsep=""), recursive = T )
  }
  cat("\n\n")
}

## Main.
# Main function, boots script.
main <- function() {
  cat("\n\n")
  ### Call functions
  ## Settings.
  # Loads defualt settings.
  settings()

  ## parseOpts.
  # Parses command line options.
  parseOpts()
  
  ### Output
  ## Show output on terminal.
  # Let the user know the program started.
  cat("Starting eQTL-mapping.R..\n\nLoading files..\n\n")
  
  # Print if the files are loaded.
  cat("Files loaded:\n",as.character(args[2][[1]]),"\n",as.character(args[3]),
      "\n\nMapping Type:\n",as.character(args[1][[1]]),"\n")
  
  ### Running the directory preperation
  ## Loading files
  # Checking existing paths.
  cat("\nCheck if directorys exist..\n\n")
  
  ## Checking directories.
  # Main directory, should be universal on every system.
  global("mainDir", "~/Documents/")
  
  # Sub directory.
  global("subDir", "R/")
  
  ### Prepare data
  ## prepare data files for mapping.
  # Call data_prep to prepare data.
  data_prep()
  
  ## Check Type.
  # Check the mapping type.
  if ( MAPPING_TYPE == 'basic') {
    
    # Create directory for the Basic results.
    createDir("Basic")
    
    # Map eQTLs.
    eQTL()
    
  } else if ( MAPPING_TYPE == 'ct' ) {
    
    # Create directory for the Cis results.
    createDir("Cis")
    
    # Create directory for the Trans results.
    createDir("Trans")
    
    # Map cis- trans-eQTLs.
    ctQTL()
    
  } else if ( MAPPING_TYPE == 's' ) {
    
    # Create directory for the Splice results.
    createDir("Splice")
    
    # Map splice QTLs.
    sQTL()
    
  } else {
    error("No mapping type accepted", "Something unexpected")
  }
}

main()
