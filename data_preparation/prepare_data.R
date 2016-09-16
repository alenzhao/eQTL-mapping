######################
# Author: E. Schutte #
######################

### Settings
## Command line arguments
# Makes command line arguments availble for this script.
args = commandArgs(trailingOnly = T)

# Test if there is at least one argument: if not, return an exception.
if ( length(args) == 0 )  {
  stop("At least two argument must be supplied!\nRscript --vanilla prepare_data.R <genotype_file> <gene_expression_file> [output]\n", call.=F)
} else if ( length(args) <= 2 ) {
  #default output
  args[3] = "out.txt"
} else {
  stop("This is for testing, I don't know what will happen")
}
### Functions
## Change_allele_to_frequencies
# Function to identify the genotype and determine wheter it is major, minor or heterozygote.
change_allele_to_frequencies <- function (x,y) {
  if ( identical(x,allele.list[[1]]$genotype) ) {
    # The genotype is major, return 0.
    return(0)
    
  } else if ( identical(allele.list[[4]]$genotype,x) ) {
    # The genotype is NA or 00, return 9.
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
mainDir <- "~/Documents"

# Sub directory.
subDir <- "/R-Output/R/eQTL-mapping"

# If the main and sub directory do not exist, create them.
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

cat("Starting prepare_data.R..\n\nLoading files..\n\n")
load(args[1])
load(args[2])
cat("Files loaded:\n",args[1],"\n",args[2],"\n\n")

## Transform the snps
# For readability change the matrix name and remove Norway data.
snps <- mergedGenotypes_clones2[-(23:27),]

# Keeps track of the SNPs that are causing noise in the data and are not representable.
snps.discarded <- c()
snps.discarded.pos <- c()
snps.discarded.counter <- 0

# For every column in the genotype data.
for (i in 1:ncol(snps)) {
  
  # Print current column and column name.
  cat("#",i,"- current column: ",colnames(snps)[i],"\n")
  
  # Save genotype and genotype count for each column.
  alleles.all <- table(snps[,i])
  
  # A minimum of 3 samples for a genotype is required to effectively calculate correlation.
  # Therefore we are removing the rows where this is not the case.
  if ( length(table(snps[,i])) == 1 | min(table(snps[,i])) < 3 ) {
    
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
    
    cat("Current allele = ", geno,'\n')
    
    # Extract the number of occurrences.
    count <- alleles.all[[a]]
    
    # Split the genotype and check if the first allele is the same as the second
    # thus checking for heterozygote or homozygote genotype.
    alleles <- strsplit(geno, "")
    
    if (identical(alleles[[1]][1],alleles[[1]][2]) & geno != '00') {
      cat('Alleles are identical, sorted as homozygote\n')
      
      if (highest_count_homozygote == 0) {
        # For the first iteration, the first genotype's occurences are the most occured
        # genotype and will be added on the first position.
        cat("Adding ",geno," with count ", count, "as major allele\n")
        highest_count_homozygote <- count
        allele.list[[1]]$genotype <- geno
        allele.list[[1]]$count <- count
        
      } else if (highest_count_homozygote < count) {
        # If the current genotype has more occurences than the lastly noted genotype,
        # it is automatically the major allele on the first position.
        cat("Adding ",geno," with count ", count, "as major allele\nPreviously ",
            allele.list[[1]]$genotype," with count ", allele.list[[1]]$count,'\n')
        
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
        cat("Adding ",geno," with count ", count, "as minor allele\n")
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
        cat('alleles are not identical, heterozygote\n')
        cat("Adding ",geno," with count ", count, "as heterozygote allele\n")
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
    cat("Changing current snip: ",snps[j,i]," to: ", snp.freq,"\n")
    
    # Save the SNP frequency at previous occupied genotype location.
    snps[j,i] <- snp.freq
  }
}

# Remove discarded snps from data set.
snps <- snps[,-snps.discarded.pos]

# Convert genotype matrix to numeric values for Matrix eQTL analaysis.
class(snps) <- "numeric"

# Transfer the genotype matrix so that the columns 'align' with the gene expression matrix.
snps.t <- t(snps)

# Remove those samples that don't occur in the gene expression data. 
# “TCC-09-1" and ”TCC-08-1”.
indices.snps <- c()
for (colname in colnames(count.all)){
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
  indices.ge <- c(indices.ge,  grep(colname, colnames(count.all)))
}

# Re-arange data frame
count.all <- count.all[,indices.ge]

### Save image
## Save the current data imgae for further use.
# Save this particular data image for basic eqtl mapping.
save.image(file=paste(getwd(),"/data_preparation/RData/basic_eqtl_mapping.RData",sep=""))
save(count.all,snps.t, file=paste(getwd(),"/data_preparation/RData/basic_eqtl_mapping.RData",sep=""))