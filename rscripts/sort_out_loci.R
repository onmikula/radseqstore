### ARGS 
# 1 - sequences, e.g. 'heterocephalus.alleles' 
# 2 - format, e.g. 'ipyrad'
# 3 - bedfile mapping the loci to reference genome
# 4-... - bedfiles defining one or more classes of loci, e.g. 'heterocephalus_coding.bed', 'heterocephalus_noncoding.bed' 

args <- commandArgs(TRUE)

### TOOLS
library(Biostrings)
library(ape)
source("read_loci.R")
source("write_loci.R")
read_loci <- match.fun(paste0("read", "_", tolower(args[2])))


### SORTING OUT LOCI
loci <- read_loci(args[1], class="DNAbin")

refbed <- read.delim(args[3], header=FALSE)
bedfiles <- args[-c(1,2,3)]

for (i in seq_along(bedfiles)) {
	bed <- read.delim(bedfiles[i], header=FALSE)
	subset <- base::intersect(refbed[,4], bed[,4])
	bed <- bed[bed[,4] %in% subset,]
	start <- refbed[match(bed[,4], refbed[,4]),2]
	bed[,2] <- bed[,2] - start + 1
	bed[,3] <- bed[,3] - start
	loc <- lapply(loci[subset], as.matrix)
	for (j in seq_along(subset)) {
		k <- unlist(lapply(which(bed[,4] == subset[j]), function(jj) bed[jj,2]:bed[jj,3]))
		loc[[j]] <- loc[[j]][,k]
	}
	nexus_file <- sub("\\.bed$", ".nexus", bedfiles[i])
	write_nexus(loc, file=nexus_file)
}


