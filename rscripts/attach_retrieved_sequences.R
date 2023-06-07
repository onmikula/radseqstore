### ARGS 
# 1 - retrievedseq, e.g. 'heterocephalus_outgroup-loci.fasta'
# 2 - sequences, e.g. 'heterocephalus.alleles' 
# 3 - format, e.g. 'ipyrad'

args <- commandArgs(TRUE)

### TOOLS
library(Biostrings)
library(ape)
source("read_loci.R")
read_loci <- match.fun(paste0("read", "_", tolower(args[3])))


### INPUT DATA
loci <- read_loci(args[2], class="DNAStringSet")
retrievedseq <- setNames(Biostrings::readDNAStringSet(args[1], format="fasta"), names(loci))

### ATTACHING SEQUENCES
loci <- lapply(loci, function(x) strsplit(as.character(x), ""))
retrievedseq <- lapply(retrievedseq, function(x) strsplit(as.character(x), ""))
for (i in seq_along(retrievedseq)) names(retrievedseq[[i]]) <- "ref"
dir.create("tobejoined")
for (s in names(loci)) {
	ape::write.dna(ape::as.DNAbin(c(retrievedseq[[s]], loci[[s]])), file=paste0("tobejoined/", s, ".fasta"), format="fasta")
}
writeLines(names(loci), "tobejoined.txt")

