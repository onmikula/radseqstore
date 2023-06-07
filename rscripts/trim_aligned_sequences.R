### ARGS 
# 1 - reference sequences, e.g. 'heterocephalus_annotated-loci.fasta' 
# 2 - bed table of reference sequences, e.g. 'heterocephalus.bed' 

args <- commandArgs(TRUE)


### TOOLS
library(Biostrings)
library(ape)


### SEQUENCES
alignedfiles <- list.files("tobejoined", full.names=TRUE)
aligned <- setNames(lapply(alignedfiles, ape::read.dna, format="fasta"), gsub("^tobejoined/|\\..+$", "", alignedfiles))
bed <- read.delim(args[2], header=FALSE)
refer <- setNames(Biostrings::readDNAStringSet(args[1], format="fasta", use.names=TRUE), bed[,4])[names(aligned)]

### TRIMMING
for (i in seq_along(aligned)) {
	ann <- which(rownames(aligned[[i]]) == "ref")
	pos <- regexpr(as.character(refer[[i]]), paste(aligned[[i]][ann,], collapse=""))
	pos <- pos - 1 + seq(attr(pos,"match.length"))
	ape::write.dna(aligned[[i]][,pos], alignedfiles[i], format="fasta")
}

