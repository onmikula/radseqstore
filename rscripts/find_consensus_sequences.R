### ARGS 
# 1 - sequences, e.g. 'heterocephalus.alleles'
# 2 - format, e.g. 'ipyrad'

args <- commandArgs(TRUE)


### TOOLS
library(Biostrings)
source("read_loci.R")
read_loci <- match.fun(paste0("read", "_", tolower(args[2])))


### INPUT
loci <- read_loci(args[1], class="DNAStringSet")


### OUTPUT
outputname <- paste0(gsub("^.*\\/|\\.[[:alpha:]]+$", "", args[1]), "_consensus.fasta")
consmatrices <- lapply(loci, function(x) Biostrings::consensusMatrix(x)[names(IUPAC_CODE_MAP),])
consensus <- lapply(consmatrices, Biostrings::consensusString, ambiguityMap=IUPAC_CODE_MAP, threshold=0.25)
writeLines(as.character(rbind(paste0(">", names(consensus)), consensus)), con=outputname)

