### ARGS 
# 1 - retrievedseq, e.g. 'heterocephalus_outgroup-loci.fasta'
# 2 - sequences, e.g. 'heterocephalus.alleles' 
# 3 - format, e.g. 'ipyrad'
# 4 - bedfile, e.g. 'heterocephalus_consensus.bed'

args <- commandArgs(TRUE)

### TOOLS
library(Biostrings)
library(ape)
source("read_loci.R")
read_loci <- match.fun(paste0("read", "_", tolower(args[3])))


### INPUT DATA
loci <- read_loci(args[2], class="DNAStringSet")


### OUTPUT OF MAPPING
retrievedseq <- Biostrings::readDNAStringSet(args[1], format="fasta")
if (!is.na(args[4])) {
	bedtable <- read.delim(args[4], header=FALSE)
	gencoords <- names(retrievedseq)
	names(retrievedseq) <- bedtable[,4]
	single <- !gencoords %in% gencoords[duplicated(gencoords)]
	retrievedseq <- retrievedseq[single]
	bedtable <- bedtable[single,]
	write.table(bedtable, file=sub("\\.bed$", "_filtered.bed", args[4]), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
} else {
	names(retrievedseq) <- names(loci)
}


### PREPARING FOR ALIGNMENT
loci <- lapply(loci[names(retrievedseq)], function(x) strsplit(as.character(x), ""))
retrievedseq <- lapply(retrievedseq, function(x) strsplit(as.character(x), ""))
for (i in seq_along(retrievedseq)) names(retrievedseq[[i]]) <- "ref"
dir.create("tobejoined")
for (s in names(retrievedseq)) {
	ape::write.dna(ape::as.DNAbin(retrievedseq[[s]]), file=paste0("tobejoined/", s, "_ref.fasta"), format="fasta")
	ape::write.dna(ape::as.DNAbin(loci[[s]]), file=paste0("tobejoined/", s, "_seq.fasta"), format="fasta")
}
writeLines(names(retrievedseq), "tobejoined.txt")

