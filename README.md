# radseqstore
A storage room for scripts and functions used for analyses of (dd)RAD sequences.

The bash scripts were written for computations at machines of MetaCentrum, the Czech academic computational grid (https://metavo.metacentrum.cz). They can contain some features, that are specific to this computational grid and, thus, they may need to be modified before being used elsewhere. The same applies to paths included in the scripts as they don't match the directory structure of this repository. All the scripts are shortly commented below.

**seq-from-reference-genome**
takes alignments of RAD loci and the reference genome that was used for their assembly and for each of the loci it extracts the corresponding part of the reference genome and attaches these sequences to the alignments. The first argument of the script is the reference genome (either a gzipped file or a directory containing genome already indexed by bwa software). The second argumnet is an output of RADseq assembly in the for of full phased sequences. The third arguments specifies the software used for the assembly and hence the format of the sequence file - 'ipyrad' or 'stacks' are implemented. The script also retrieves genomic position of the loci (in oordinates of the reference genome) and outputs them in the form bed table.
