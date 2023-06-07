### ARGUMENTS 
# 1 - output bedfile with shared loci
# 2 - file with names of input bedfiles
# 3 - file with names of subset bedfiles
# 4 - file with names of nexus files

args <- commandArgs(TRUE)
outputbed <- args[1]
inputbeds <- readLines(args[2])
subsetbeds <- readLines(args[3])
nexfiles <- readLines(args[4])


### TOOLS
library(ape)
source("read_loci.R")
source("write_loci.R")

shared <- read.delim(outputbed, header=FALSE)
shared <- shared[!duplicated(shared[,c(1:3,6)]),]
shared[,7] <- paste(shared[,1], paste(shared[,2], shared[,3], sep="-"), shared[,6], sep=":")

no <- length(inputbeds)
loci <- setNames(vector("list", no), gsub("^.*\\/|\\.bed$", "", inputbeds))
for (i in seq(no)) {
	subset <- read.delim(subsetbeds[i], header=FALSE)
	subset[,7] <- paste(subset[,1], paste(subset[,2], subset[,3], sep="-"), subset[,6], sep=":")
	subset <- subset[match(shared[,7], subset[,7]),]
	input <- read.delim(inputbeds[i], header=FALSE)
	input <- input[match(subset[,4], input[,4]),]
	loci[[i]] <- read_nexus(nexfiles[i], class="matrix")[subset[,4]]
	for (j in seq_along(loci[[i]])) {
		jj <- (subset[j,2] - input[j,2] + 1):(subset[j,3] - input[j,2])
		loci[[i]][[j]] <- loci[[i]][[j]][,jj]
	}
}

combined <- setNames(vector("list", nrow(shared)), shared[,7])
for (j in seq_along(combined)) {
	combined[[j]] <- do.call(rbind, lapply(loci, "[[", j))
}

write.table(shared, outputbed, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write_nexus(combined, sub("\\.bed$", ".nexus", outputbed))

