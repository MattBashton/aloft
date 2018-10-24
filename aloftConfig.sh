#!/bin/bash -eu

############################
# Setting for aloft script #
############################

# Read length for STAR index, default here is 201bp
# NOTE LEAVE AS IS SEE HERE:
# https://github.com/suhrig/arriba/issues/9
READL=201
# Which = 200bp for use in actual script

# CPU threads for GNU parallel
GNUP_THREAD=8
# CPU threads for STAR index build
STARI_CPU=8
# CPU thread for STAR alignment
STARA_CPU=8
# CPU thread for GNU Parralle concurrent STAR alignments, note above setting!
# Default here = 2 x 8 thread STAR alignment jobs = 16 threds in total
# Note each STAR job = ~ 30GB of RAM
GNUP_THREAD_STAR=2

# Arriba release
ARRIBA_VER="v1.0.1"

# Inferred from above
ARRIBA_DIR="arriba_${ARRIBA_VER}"
ARRIBA="${ARRIBA_DIR}/arriba"
ARRIBA_REL="https://github.com/suhrig/arriba/releases/download/${ARRIBA_VER}/arriba_${ARRIBA_VER}.tar.gz"

# Max reads for Arriba to use, it will downsample above this count
MAX_READS="300"
# Other extra Arrib arguments, place these into this string,
# left empty by default
ARRIBA_ARG=""

# Index settings
# Referance and annotaion to use
REF="GRCh38"
ANNO="ENSEMBL93"
REF_ANNO="${REF}+${ANNO}"
STAR_INDEX_DIR="STAR_index_${REF}_${ANNO}"
# Compression level for temp BAM file from STAR
# Set to zero to speed up performance
STAR_BAM_CMPLVL=0
BLACKLIST="${ARRIBA_DIR}/database/blacklist_hg38_GRCh38_2018-04-04.tsv.gz"
CYTOBANDS="${ARRIBA_DIR}/database/cytobands_hg38_GRCh38_2018-02-23.tsv"
DOMAINS="${ARRIBA_DIR}/database/protein_domains_hg38_GRCh38_2018-03-06.gff3"

# COSMIC "Complete Fusion Export" download location:
# https://cancer.sanger.ac.uk/cosmic/download
# Log-in required, please download and place in current dir
COSMIC="CosmicFusionExport.tsv.gz"

# STAR alignment options
# Cribbed from https://arriba.readthedocs.io/en/v1.0.0/workflow/
STAR_ARG="--genomeLoad NoSharedMemory --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression ${STAR_BAM_CMPLVL} --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3"

# SAMtools RAM usage
# Note that GNUP_THREAD will be used so x this by that value to work out max
# possible RAM usage current set-up will use 8x8192MB = 64GB
# Note actual usage will be smaller if input file is less than 8GB
SAMTOOLS_RAM="8192M"
# Number of threads for samtools sort, note number of samtools jobs defaults to 8
SAMTOOLS_CPU=2
SAMTOOLS_CMPLVL=6

# ungzip command for fastq.gz files
# Use pigz here otherwise use zcat
#UNGZIP="pigz -dc"
UNGZIP="zcat"
