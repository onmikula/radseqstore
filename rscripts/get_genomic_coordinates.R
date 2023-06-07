### ARGS 
# 1 - sequences / .bed file, e.g. 'heterocephalus.alleles'
# 2 - format, e.g. 'ipyrad'

args <- commandArgs(TRUE)


### GENOMIC COORDINATES
if (grepl("\\.bed$", args[1])) {

	gencoord <- read.delim(args[1], header=FALSE)
	gencoord_vector <- setNames(paste0(gencoord[,1], ":", gencoord[,2], "-", gencoord[,3]), gencoord[,4])
	regionfile <- paste0(gsub("^.*\\/|\\.[[:alpha:]]+$", "", args[1]), "_genomic-regions.txt")
	writeLines(gencoord_vector, regionfile)

} else {

	library(Biostrings)
	source("read_loci.R")
	read_loci <- match.fun(paste0("read", "_", tolower(args[2])))

	loci <- read_loci(args[1], class="DNAStringSet")
	gencoord <- attr(loci, "gencoord")
	gencoord$Start <- gencoord$Start - 1
	bedfile <- paste0(gsub("^.*\\/|\\.[[:alpha:]]+$", "", args[1]), ".bed")
	write.table(gencoord, bedfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

	gencoord_vector <- setNames(paste0(gencoord$Scaffold, ":", gencoord$Start + 1, "-", gencoord$End), gencoord$Locus)
	regionfile <- paste0(gsub("^.*\\/|\\.[[:alpha:]]+$", "", args[1]), "_genomic-regions.txt")
	writeLines(gencoord_vector, regionfile)

}
