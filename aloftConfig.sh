#!/bin/bash -eu
############################
# Setting for aloft script #
############################
# Matt Bashton 2018

# Read length
#------------
# Read length for STAR index, default here is 201bp
# NOTE LEAVE AS IS SEE HERE:
# https://github.com/suhrig/arriba/issues/9
READL=201
# Which = 200bp for use in actual script


# GNU parallel settings
# ---------------------
# CPU threads for GNU parallel, this is used for Arriba and plotting jobs.
# This value controls the number of concurrent jobs.
GNUP_THREAD=8

# CPU threads for GNU parallel concurrent STAR alignments,
# Default here = 2x 8-thread STAR alignment jobs = 16 threads in total
# Note each STAR job = ~30GB of RAM so total ~60GB.
GNUP_THREAD_STAR=2

# CPU threads for GNU Parallel concurrent SAMtools jobs. Default here
# = 4x 2-thread SAMtools jobs = 8 threads in total
# Note each SAMtools sort worker thread = ~8GB of RAM so this would use 64GB of
# RAM.
GNUP_THREAD_SAMTOOLS=4


# Arriba settings
#----------------
# Arriba release
ARRIBA_VER="v1.0.1"

# Inferred from above
ARRIBA_DIR="arriba_${ARRIBA_VER}"
ARRIBA="${ARRIBA_DIR}/arriba"
ARRIBA_REL="https://github.com/suhrig/arriba/releases/download/${ARRIBA_VER}/arriba_${ARRIBA_VER}.tar.gz"

# Max reads for Arriba to use, it will downsample above this count
MAX_READS="300"
# Other extra Arriba arguments, place these into this string,
# left empty by default
ARRIBA_ARG=""
# Datasets form Arriba database dir, these need to match STAR assembly below
BLACKLIST="${ARRIBA_DIR}/database/blacklist_hg38_GRCh38_2018-04-04.tsv.gz"
CYTOBANDS="${ARRIBA_DIR}/database/cytobands_hg38_GRCh38_2018-02-23.tsv"
DOMAINS="${ARRIBA_DIR}/database/protein_domains_hg38_GRCh38_2018-03-06.gff3"


# STAR settings
# -------------
# STAR CPU threads for alignment
STARA_CPU=8
# Reference and annotation to use for STAR index
REF="GRCh38"
ANNO="ENSEMBL93"
REF_ANNO="${REF}+${ANNO}"
STAR_INDEX_DIR="STAR_index_${REF}_${ANNO}"
# Compression level for temp BAM file from STAR
# Set to zero to speed up performance
STAR_BAM_CMPLVL=0
# STAR alignment options
# Cribbed from Arrib 1.0.1 workflow:
# https://arriba.readthedocs.io/en/v1.0.1/workflow/
STAR_ARG="--genomeLoad NoSharedMemory --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression ${STAR_BAM_CMPLVL} --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3"
# Ungzip command for fastq.gz files for reading FASTQ into STAR
# Use pigz here if you have it otherwise use zcat
#UNGZIP="pigz -dc"
UNGZIP="zcat"


# COSMIC known fusions list generation
# ------------------------------------
# COSMIC "Complete Fusion Export" download location:
# https://cancer.sanger.ac.uk/cosmic/download
# Log-in required, please download and place in current dir
COSMIC="CosmicFusionExport.tsv.gz"


# SAMtools settings
#------------------
# SAMtools RAM usage
# Note that GNUP_THREAD will be used so x this by that value to work out max
# possible RAM usage current set-up will use 4x jobs 2x threads x 8192MB per
# thread = 64GB Note actual usage will be smaller if input file is less than
# 8GB.
SAMTOOLS_RAM="8192M"
# Number of threads for SAMtools sort, note number of SAMtools jobs defaults to 4
SAMTOOLS_CPU=2
SAMTOOLS_CMPLVL=6
