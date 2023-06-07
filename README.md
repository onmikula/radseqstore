# radseqstore
A storage room for scripts and functions used for analyses of (dd)RAD sequences.

The bash scripts were written for computations at machines of MetaCentrum, the Czech academic computational grid (https://metavo.metacentrum.cz). They can contain some features, that are specific to this computational grid and, thus, they need to be modified appropriately before being used elsewhere. The same applies to paths included in the scripts as they don't match the directory structure of this repository. All the scripts are shortly commented below. The embedded R scripts may require specific packages or custom functions, but everything should be appearent from them (see folder 'rscripts'). The custom functions are provided at the same place.

**seq-from-reference-genome**
takes alignments of RAD loci and the reference genome that was used for their assembly and for each of the loci it extracts the corresponding part of the reference genome and attaches these sequences to the alignments. The argguments are: (1) the reference genome (either a gzipped file or a directory containing genome already indexed by bwa software); (2) a sequence file, i.e., an output of RADseq assembly in the for of full phased sequences; (3) software used for the assembly and hence the format of the sequence file - 'ipyrad' or 'stacks' are implemented. The script also retrieves genomic positions of the loci (in oordinates of the reference genome) as recorded in the sequence file and outputs them in the form bed table.

**seq-from-outgroup-genome**
takes alignments of RAD loci and an outgroup genome and for each of the loci it finds consensus sequence, map it to the references genome, retrieves the corresponding sequence and aligns it to its alignment. The argguments are: (1) ...
