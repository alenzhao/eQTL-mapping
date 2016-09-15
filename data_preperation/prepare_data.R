### Settings
## Command line arguments
# Makes command line arguments availble for this script.
args = commandArgs(trailingOnly = T)

# Test if there is at least one argument: if not, return an exception.
if ( length(args) == 0 )  {
  stop("At least two argument must be supplied!\nRscript --vanilla prepare_data.R genotype_file gene_expression_file [output]\n(input file).n", call.=F)
} else if ( length(args) <= 2 ) {
  #default output
  args[3] = "out.txt"
} else {
  stop("This is for testing, I don't know what will happen")
}
