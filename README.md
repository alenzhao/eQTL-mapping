# eQTL-mapping
Software to map eQTLs, (cis- and trans-eQTLs), as well as sQTLs by comparing gene expression data with genotype data with R's Matrix eQTL package. This software, as a whole, can:
- map basic eQTLs
- map cis- and trans-eQTLs
- map sQTLs
- prepare data sets for Matrix eQTL input
- calculate MAF

# Set-up
Your gene expression data file and your genotype data file should be prepared for input into the Matrix eQTL package. By doing so, the data created will be saved in ‘/data_preparation/RData/‘.

How do we prepare data?
* Open a terminal in the R project (root folder).
* Run the following command: ‘Rscript --vanilla ./data_preparation/prepare_data.R “<your_genotype_data>“ “<your_gene_expression_data>”’ 
* The resulting *.RData file will appear in the ‘/RData/‘ folder.
* Next your SNP and Gene position files should be set-up. They have to have the exact format described here for SNP positions: ‘http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/snpsloc.txt' and here for Gene positions: ‘http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/geneloc.txt'

# Mapping basic eQTLs
* Run: ‘Rscript --vanilla ./analysis/basic_eQTL.R’
* The results are now in your ‘/Users/<your_name>/Documents/R-Output/basic’ folder.

# Mapping cis- and trans-eQTLs
* Run: ‘Rscript --vanilla ./analysis/cis-trans-eQTL.R "~/Dropbox/Erik/genepos.txt" "~/Documents/CeD_43loci_alternative_format.txt”’
* The results are now in your ‘‘/Users/<your_name>/Documents/R-Output/cis’ & ‘/Users/<your_name>/Documents/R-Output/trans’ folder.
