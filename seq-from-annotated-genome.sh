#!/bin/bash
#PBS -l walltime=96:00:00
#PBS -l select=1:ncpus=4:mem=96gb:scratch_local24gb
#PBS -N seq-from-annotated-genome

# loading modules
module add r/4.1.3-gcc-10.2.1-6xt26dl
module add bwa-0.7.17
module add samtools-1.10
module add bedtools-2.26.0
module add mafft/7.487

R_LIBS="/storage/brno7-cerit/projects/fishery/software/R/libs:$R_LIBS"
export R_LIBS

# loading R scripts & functions
rscripts_dir="/storage/brno7-cerit/projects/fishery/scripts/job/onmikula/molerats/rscripts"

find_consensus_sequences="$rscripts_dir/find_consensus_sequences.R"
cp "$find_consensus_sequences" "$SCRATCHDIR"
find_consensus_sequences="${find_consensus_sequences/*\//}"

get_genomic_coordinates="$rscripts_dir/get_genomic_coordinates.R"
cp "$get_genomic_coordinates" "$SCRATCHDIR"
get_genomic_coordinates="${get_genomic_coordinates/*\//}"

prepare_retrieved_sequences="$rscripts_dir/prepare_retrieved_sequences.R"
cp "$prepare_retrieved_sequences" "$SCRATCHDIR"
prepare_retrieved_sequences="${prepare_retrieved_sequences/*\//}"

trim_aligned_sequences="$rscripts_dir/trim_aligned_sequences.R"
cp "$trim_aligned_sequences" "$SCRATCHDIR"
trim_aligned_sequences="${trim_aligned_sequences/*\//}"

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
   # copy indexed outgroup genome to scratch
   cp "$genome"/*.* "$SCRATCHDIR"
   genomeid="${genome/*\//}"
   genomefile="$genomeid".fna
else
   # copy zipped genome to scratch, unzip it and index it
   cp "$genome" "$SCRATCHDIR"
   genome="${genome/*\//}"
   gunzip "$SCRATCHDIR/$genome"
   genomefile="${genome/.gz/}"
   bwa index "$genome"
   genomeid="${genomefile/.fna/}"
fi


# copy sequences of genomic (e.g. ddrad) loci to scratch
cp "$seqfile" "$SCRATCHDIR"
seqfile="${seqfile/*\//}"

# go into the "$SCRATCHDIR"
cd "$SCRATCHDIR"

# find or read in consenus sequences
if [ -z "$4" ]
then
   # find consensus sequences of ddrad loci
   Rscript --vanilla "$find_consensus_sequences" "$seqfile" "$format"
   conseqfile="${seqfile/.*/_consensus.fasta}"
else
   conseqfile="$4"
   cp "$conseqfile" "$SCRATCHDIR"
   conseqfile="${conseqfile/*\//}"
fi

# file names
bamfile="${seqfile/.*/.bam}"
unsortedbam="${bamfile/.bam/-unsorted.bam}"
sortedbam="${bamfile/.bam/-sorted.bam}"
bedfile="${seqfile/.*/.bed}"
annotatedseq="${seqfile/.*/_annotated-loci.fasta}"
filteredbed="${seqfile/.*/_filtered.bed}"
regionfile="${filteredbed/.*/_genomic-regions.txt}"
nexfile="${seqfile/.*/.nexus}"


# map ddrad consensus sequences to the indexed outgroup genome
# 2304: samtools flags SECONDARY,SUPPLEMENTARY (primary = neither 0x100 nor 0x800 bit set)
bwa mem "$genomefile" "$conseqfile" -t 4 | samtools view -b -o "$unsortedbam"
samtools sort -m 3500M -@ 2 -T tmp -O bam -o "$sortedbam" "$unsortedbam"
samtools view -F 2304 -b -o "$bamfile" "$sortedbam"

# get outgroup sequences: .bam -> .bed -> .fasta
# (-s option means reverse complementing of sequences from the negative strand)
bedtools bamtobed -i "$bamfile" > "$bedfile"
bedtools getfasta -s -fi "$genomefile" -bed "$bedfile" -fo "$annotatedseq"

# attach retrieved sequences
Rscript --vanilla "$prepare_retrieved_sequences" "$annotatedseq" "$seqfile" "$format" "$bedfile"

# get genomic coordinates of ddrad loci
Rscript --vanilla "$get_genomic_coordinates" "$filteredbed"

# align annotated sequences by MAFFT
while read loc
do
mafft --add "tobejoined/$loc""_ref.fasta" "tobejoined/$loc""_seq.fasta" > "tobejoined/$loc"".fasta"
rm "tobejoined/$loc""_ref.fasta"
rm "tobejoined/$loc""_seq.fasta"
done < tobejoined.txt

# trim aligned sequences
Rscript --vanilla "$trim_aligned_sequences" "$annotatedseq" "$bedfile"

# join aligned sequences
Rscript --vanilla "$join_aligned_sequences" "$seqfile" "$format" "$genomeid"

# move output files to output directory
if [ ! -d output ]
then
  mkdir output
fi

mv "$conseqfile" "output/$conseqfile"
mv "$bamfile" "output/$bamfile"
mv "$bedfile" "output/$bedfile"
mv "$filteredbed" "output/$filteredbed"
mv "$regionfile" "output/$regionfile"
mv "$annotatedseq" "output/$annotatedseq"
mv "$nexfile" "output/$nexfile"
