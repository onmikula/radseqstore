### FUNCTION
# write_nexus
### ARGUMENTS
# loci - input file with phylogenetic trees or DNA alignment
#	trees: object of class "phylo", "multiPhylo" or a list of "phylo" objects
#	DNA: character matrix or object of class "DNAbin" 
# file - name of output file

write_nexus <- function(loci, file, part=NULL, missing="-", gap="!") {

	if (is.null(part)) {
		part <- attr(loci, "part")
	}
	if (!is.list(loci) | inherits(loci, "DNAbin")) {
		loci <- setNames(list(loci), "sequence")
	}
	if (inherits(loci[[1]], "matrix")) {
		nbp <- sapply(loci, ncol)
		seqnames <- sort(unique(unlist(lapply(loci, rownames))))
	} else if (inherits(loci[[1]], "DNAbin")) {
		nbp <- sapply(loci, ncol)
		seqnames <- sort(unique(unlist(lapply(loci, rownames))))
	} else if (inherits(loci[[1]], "DNAStringSet")) {
		nbp <- sapply(loci, function(x) nchar(x[1]))
		seqnames <- sort(unique(unlist(lapply(loci, names))))		
	}
	if (is.null(part)) {
		cs <- cumsum(nbp)
		part <- matrix(c(1, cs[-length(nbp)] + 1, cs), length(nbp), 2, dimnames=list(names(nbp), NULL))
	}	
	if (inherits(loci[[1]], "matrix")) {
		empty <- lapply(nbp, function(x) rep(missing, x))
	} else if (inherits(loci[[1]], "DNAbin")) {
		empty <- lapply(nbp, function(x) ape::as.DNAbin(rep(missing, x)))
	} else if (inherits(loci[[1]], "DNAStringSet")) {
		empty <- lapply(nbp, function(x) Biostrings::DNAStringSet(paste(rep(missing, x), collapse="")))
	}

	sequences <- setNames(vector("list", length(seqnames)), seqnames)
	if (inherits(loci[[1]], "matrix")) {
		for (i in seq_along(seqnames)) {
			matches <- sapply(lapply(loci, rownames), function(x) match(seqnames[i], x))
			sequences[[i]] <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][matches[j],] else character(0))
			loci <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][-matches[j],,drop=FALSE] else loci[[j]])
			zero <- sapply(sequences[[i]], length) == 0
			sequences[[i]][zero] <- empty[zero]
			sequences[[i]] <- paste(unlist(sequences[[i]]), collapse="")
		}		
	} else if (inherits(loci[[1]], "DNAbin")) {
		for (i in seq_along(seqnames)) {
			matches <- sapply(lapply(loci, rownames), function(x) match(seqnames[i], x))
			sequences[[i]] <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][matches[j],] else character(0))
			loci <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][-matches[j],,drop=FALSE] else loci[[j]])
			zero <- sapply(sequences[[i]], length) == 0
			sequences[[i]][zero] <- empty[zero]
			sequences[[i]] <- paste(unlist(lapply(lapply(sequences[[i]], as.character), as.vector)), collapse="")
		}		
	} else if (inherits(loci[[1]], "DNAStringSet")) {
		for (i in seq_along(seqnames)) {
			matches <- sapply(lapply(loci, names), function(x) match(seqnames[i], x))
			sequences[[i]] <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][matches[j],] else character(0))
			loci <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][-matches[j],,drop=FALSE] else loci[[j]])
			zero <- sapply(sequences[[i]], length) == 0
			sequences[[i]][zero] <- empty[zero]
			sequences[[i]] <- paste(unlist(lapply(sequences[[i]], as.character)), collapse="")
		}		
	}

	ntaxa <- length(seqnames)
	nsites <- sum(nbp)
	format <- paste("FORMAT", "DATATYPE=DNA", paste0("GAP=", gap), paste0("MISSING=", missing, ";"))
	header <- c("#NEXUS\n", "BEGIN DATA;", paste0("DIMENSIONS", "  ", "NTAX=", ntaxa, " ", "NCHAR=", nsites, ";"), format, "MATRIX\n")
	sequences <- unname(c(seqnames, unlist(sequences))[rep(seq(ntaxa), each=2) + (seq(2 * ntaxa) %% 2 == 0) * ntaxa])			
	nexus <- c(header, sequences, ";\n", "END;")
	if (!is.null(part)) {
		p <- apply(part, 1, function(x) paste(x, collapse=" - "))
		p <- paste(names(p), p, sep=" = ")
		p <- paste(paste(paste("\t", "CHARSET", sep=""), p), ";", sep="")
		p <- c("BEGIN assumptions;", p, "END;")
		nexus <- c(nexus, "\n", p)
	}
	
	writeLines(nexus, con=file, sep="\n", useBytes=FALSE)

}

