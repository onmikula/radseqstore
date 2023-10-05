#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=4:mem=96gb:scratch_local=12gb
#PBS -N seq-from-reference-genome

# loading modules
module add r/4.1.3-gcc-10.2.1-6xt26dl
module add samtools-1.10
module add mafft/7.487

R_LIBS="/storage/brno7-cerit/projects/fishery/software/R/libs:$R_LIBS"
export R_LIBS

# loading R scripts & functions
rscripts_dir="/storage/brno7-cerit/projects/fishery/scripts/job/onmikula/molerats/rscripts"

get_genomic_coordinates="$rscripts_dir/get_genomic_coordinates.R"
cp "$get_genomic_coordinates" "$SCRATCHDIR"
get_genomic_coordinates="${get_genomic_coordinates/*\//}"

attach_retrieved_sequences="$rscripts_dir/attach_retrieved_sequences.R"
cp "$attach_retrieved_sequences" "$SCRATCHDIR"
attach_retrieved_sequences="${attach_retrieved_sequences/*\//}"

join_aligned_sequences="$rscripts_dir/join_aligned_sequences.R"
cp "$join_aligned_sequences" "$SCRATCHDIR"
join_aligned_sequences="${join_aligned_sequences/*\//}"

read_loci="$rscripts_dir/read_loci.R"
cp "$read_loci" "$SCRATCHDIR"

write_loci="$rscripts_dir/write_loci.R"
cp "$write_loci" "$SCRATCHDIR"

# arguments
genome="$1"
seqfile="$2"
format="$3"

if [ -d "$genome" ]
then
   genomefile="$genome/${genome/*\//}.fna"
   cp "$genomefile" "$SCRATCHDIR"
   genomefile="${genomefile/*\//}"
else
   # copy zipped genome to scratch and unzip it
   cp "$genome" "$SCRATCHDIR"
   genome="${genome/*\//}"
   gunzip "$SCRATCHDIR/$genome"
   genomefile="${genome/.gz/}"
fi
genomeid="${genomefile/.fna/}"


# copy sequences of genomic (e.g. ddrad) loci to scratch
cp "$seqfile" "$SCRATCHDIR"

# file names
seqfile="${seqfile/*\//}"
bedfile="${seqfile/.*/}".bed
regionfile="${seqfile/.*/}"_genomic-regions.txt
referenceseq="${seqfile/.*/}"_reference-loci.fasta
alignedseq="${seqfile/.*/.nexus}"

# go into the "$SCRATCHDIR"
cd "$SCRATCHDIR"

# get genomic coordinates of ddrad loci
Rscript --vanilla "$get_genomic_coordinates" "$seqfile" "$format"

# get reference sequences corresponding to the coordinates
samtools faidx --mark-strand sign "$genomefile" -r "$regionfile" -o "$referenceseq"

# attach retrieved sequences
Rscript --vanilla "$attach_retrieved_sequences" "$referenceseq" "$seqfile" "$format"

# join aligned sequences
Rscript --vanilla "$join_aligned_sequences" "$seqfile" "$format" "$genomeid"

# move output files to output directory
if [ ! -d output ]
then
  mkdir output
fi

mv "$bedfile" "output/$bedfile"
mv "$regionfile" "output/$regionfile"
mv "$referenceseq" "output/$referenceseq"
mv "$alignedseq" "output/$alignedseq"
