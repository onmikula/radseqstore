### ARGS 
# 1 - sequences, e.g. 'heterocephalus.alleles' 
# 2 - format, e.g. 'ipyrad'
# 3 - genome, e.g. "Fukomys_damarensis_GCA_012274545.1_DMR_v1.0_HiC_genomic"

args <- commandArgs(TRUE)


### TOOLS
library(Biostrings)
library(ape)
source("read_loci.R")
source("write_loci.R")
read_loci <- match.fun(paste0("read", "_", tolower(args[2])))


### JOINING LOCI
loci <- read_loci(args[1], class="matrix")
alignedfiles <- list.files("tobejoined", full.names=TRUE)
aligned <- setNames(lapply(alignedfiles, ape::read.dna, format="fasta"), gsub("^tobejoined/|\\..+$", "", alignedfiles))
mapped <- match(names(aligned), names(loci))
loci[mapped] <- aligned
for (i in mapped) {
	r <- which(rownames(loci[[i]]) == "ref")
	rownames(loci[[i]])[r] <- args[3]
}
nexus_file <- sub("\\.[[:alpha:]]+$", ".nexus", args[1])
write_nexus(loci, file=nexus_file)
filtered_nexus_file <- sub("\\.[[:alpha:]]+$", "_filtered.nexus", args[1])
write_nexus(loci[names(loci) %in% names(aligned)], file=filtered_nexus_file)
