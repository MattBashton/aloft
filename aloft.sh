#!/bin/bash -eu
tput bold
echo "Aloft a script for running Arriba RNA-Seq fusion detection"
echo "Matt Bashton 2018"
echo " "
tput sgr0

# We need a sample sheet to run
[ $# -ne 3 ] && { echo -en " *** Error Nothing to do, usage: $(basename $0) <tab delimited sample sheet> <input FASTQ path> <output dir for run> ***\n\n" ; exit 1; }

# Set up bash to fail if pipe fails
set -o pipefail

# Where are we and what is the time?
hostname
date
BASE_DIR="$PWD"
echo -ne "\n"

# Source our config file
source aloftConfig.sh

# Work out number of samples in sample_sheet.txt
SAMPLE_SHEET="$1"
INPUT_DIR="$2"
OUTPUT_DIR="$3"
COUNTER=1
SAMPLES=$(wc -l ${SAMPLE_SHEET} | awk '{print $1}')

# Output config
echo -ne "Config\n"
echo " * Working dir: ${BASE_DIR}"
echo " * Input FASTQ dir: ${INPUT_DIR}"
echo " * Read Lenght: ${READL}bp (please ignore this - for testing only)"
echo " * Sample count: $SAMPLES"
echo " * Arriba version: ${ARRIBA_VER}"
echo " * Arriba release: ${ARRIBA_REL}"
echo " * Arriba dir: ${ARRIBA_DIR}"
echo " * Arriba binary: ${ARRIBA}"
echo " * Genome ref to use: ${REF}"
echo " * Gene annotaion to use: ${ANNO}"
echo " * Max reads for Arriba to use: ${MAX_READS}"
echo " * Optional Arriba arguments: ${ARRIBA_ARG}"
echo " * STAR index dir: ${STAR_INDEX_DIR}"
echo " * STAR ungzip command: ${UNGZIP}"
echo " * STAR BAM output (temp) compression level: ${STAR_BAM_CMPLVL}"
echo " * Concurrent STAR jobs: ${GNUP_THREAD_STAR}"
echo " * STAR alignment threads ${STARA_CPU}"
echo " * Blacklist ot use: ${BLACKLIST}"
echo " * Cytobands to use for plots: ${CYTOBANDS}"
echo " * Pfam domains to use for plots: ${DOMAINS}"
echo " * COMSIC Complete fusion report file: ${COSMIC}"
echo " * Final output dir set to: ${OUTPUT_DIR}"
echo " * Concurrent SAMtools jobs: ${GNUP_THREAD_SAMTOOLS}"
echo " * RAM allocated for SAMtools sort: ${SAMTOOLS_RAM}"
echo " * CPU threads for SAMTtools sort: ${SAMTOOLS_CPU}"
echo " * SAMtools compression level for BAM: ${SAMTOOLS_CMPLVL}"
echo -ne "\n"

# Detect if Arriba is already downloaded and built
echo "Detecting if Arriba is present"
if [[ -f ${ARRIBA} ]]
then
    echo -en " * Arriba binary found, no need to fetch and compile\n"
else
    echo -en " * Arriba not found, downloading and compiling from source...\n"
    # rm incomplete downloads or extractions
    [[ -e arriba_${ARRIBA_VER}.tar.gz ]] && rm arriba_${ARRIBA_VER}.tar.gz
    [[ -e arriba_${ARRIBA_VER} ]] && rm -rf arriba_${ARRIBA_VER}
    wget -q --show-progress ${ARRIBA_REL}
    tar -xf arriba_${ARRIBA_VER}.tar.gz
    cd ${ARRIBA_DIR}
    make
    cd $BASE_DIR
fi

# Detect if STAR index + referance and annotation is present
echo -ne "\nDetecting if STAR index is present\n"
if [[ -d ${STAR_INDEX_DIR} ]]
then
    echo -en " * STAR index ${STAR_INDEX_DIR} found, no need to rebuild\n"
else
    echo -en " * STAR index not found, downloading and building on $STARI_CPU cpu core(s)... (this may take some time)\n"
    echo -en " * Adjusting build script for a read length of ${READL}bp\n"
    
    # Edit download script with sed
    LENGTH=$(($READL-1))
    sed "s/sjdbOverhang 200/sjdbOverhang ${LENGTH}/" ${ARRIBA_DIR}/download_references.sh > ${ARRIBA_DIR}/dl_ref.sh
    chmod +x ${ARRIBA_DIR}/dl_ref.sh
    # Run script with options
    ${ARRIBA_DIR}/dl_ref.sh "${REF}+${ANNO}" ${STARI_CPU}
fi

# Remove output dir file if exists
if [[ -d ${OUTPUT_DIR} ]]
then
    rm -rf ${OUTPUT_DIR}
fi

# Make output dir 
mkdir ${OUTPUT_DIR}

# Remove Logs dir if file exists
# STAR Logs will be placed here
if [[ -d logs ]]
then
    rm -rf logs
fi

# Make Logs dir
mkdir logs

# Now Generate STAR alignment jobs
echo -ne "\nGenerating STAR alignment jobs...\n"
echo -ne " * Threads to use for STAR alignment: $STARA_CPU\n"
echo -ne " * STAR index dir is: $STAR_INDEX_DIR\n"

# Remove existing jobs file if exists
if [[ -e STAR_jobs.txt ]]
then
    rm STAR_jobs.txt
fi

# Loop over samples create STAR jobs
while [[ $COUNTER -le $SAMPLES ]];
do    
    echo -ne "  - Working on sample $COUNTER of $SAMPLES, ID: "
    # Get input files for this run
    R1=$(cut -f2- ${SAMPLE_SHEET} | awk "NR==${COUNTER}" | tr '\t' '\n' | grep R1 | awk -v var="${INPUT_DIR}" '{printf(var"/%s,",$0)}' | sed 's/,\s*$//')
    R2=$(cut -f2- ${SAMPLE_SHEET} | awk "NR==${COUNTER}" | tr '\t' '\n' | grep R2 | awk -v var="${INPUT_DIR}" '{printf(var"/%s,",$0)}' | sed 's/,\s*$//')
    # echo "$R1 $R2"
    # Set output BAM name, we want to save this to view in IGV later so saving to disk rather than pipe
    SAMPLE_ID=$(awk "NR==${COUNTER}" ${SAMPLE_SHEET} | cut -f 1)
    echo -ne "${SAMPLE_ID}"
    OUTPUT="${SAMPLE_ID}.bam"
    echo -ne ", STAR output will be: $SAMPLE_ID.bam"
    # Rest of command-line settings are taken from aloftConfig.sh
    echo "STAR --runThreadN ${STARA_CPU} --genomeDir ${STAR_INDEX_DIR} --readFilesCommand ${UNGZIP} --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE_ID} --outTmpDir TmpDir_${SAMPLE_ID} --readFilesIn ${R1} ${R2} ${STAR_ARG} > ${OUTPUT_DIR}/${SAMPLE_ID}.bam" >> STAR_jobs.txt
    echo -ne "\n"
    ((COUNTER++))
done

# Generate COSMIC known fusions list
if [[ -e CosmicFusionsList.tsv ]]
then
    COSMIC_FUSIONS=$(wc -l CosmicFusionsList.tsv | cut -d ' ' -f 1)
    echo -ne "\n${COSMIC_FUSIONS} known COSMIC fusions found in CosmicFusionsList.tsv, skipping file regeneration\n"
else
    echo -ne "\nGenerating COSMIC known fusions list\n"
    echo -ne " * Parsing: ${COSMIC} and reprocessing to CosmicFusionsList.tsv...\n"
    zcat ${COSMIC} | grep -v 'Sample ID' | cut -f 12 | grep -v -e '^$' | perl -lne '/^(\S+?)\{\S+\d+}\S+_(\S+?)\{/; print "$1\t$2"' | sort -u > TempCosmicFusionsList.tsv
    COSMIC_FUSIONS=$(wc -l TempCosmicFusionsList.tsv | cut -d ' ' -f 1)
    echo -ne " * Found ${COSMIC_FUSIONS} known events\n"
    echo -ne " * Fixing gene name issues in CosmicFusionsList.tsv...\n"
    # Remove one bogus pair, and remove temp file
    grep -v '2229+2522ALK' TempCosmicFusionsList.tsv > CosmicFusionsList.tsv
    rm TempCosmicFusionsList.tsv
    # String for sed
    # format s/OLD/NEW/g;
    # Make array 
    declare -A SED_ARRAY
    # Populate key value pairs
    SED_ARRAY[C15orf55]=NUTM1
    SED_ARRAY[C2orf44]=WDCP
    SED_ARRAY[CTAGE5]=MIA2
    SED_ARRAY[ERO1L]=ERO1A
    SED_ARRAY[LHFP]=LHFPL6
    SED_ARRAY[HN1]=JPT1
    SED_ARRAY[KIAA1598]=SHTN1
    SED_ARRAY[CASC5]=KNL1
    SED_ARRAY[KIAA0284]=CEP170B
    SED_ARRAY[KIAA1524]=CIP2A
    SED_ARRAY[MLLT4]=AFDN
    SED_ARRAY[FAM22A]=NUTM2A
    SED_ARRAY[FAM22B]=NUTM2B
    # Generate string for sed in loop
    SED_STRING=""
    for OLD in "${!SED_ARRAY[@]}"
    do
	echo -ne "  - Replace ${OLD} with: "
	echo -ne "${SED_ARRAY[${OLD}]}\n"
	STRING_TO_ADD="s/${OLD}/${SED_ARRAY[${OLD}]}/g; "
       	SED_STRING=${SED_STRING}${STRING_TO_ADD}
    done
    # Remove final two charaters
    SED_STRING=${SED_STRING:0:-2}
    # Add closing single quote
    #SED_STRING="${SED_STRING}'"
    echo -ne " * Final string is: ${SED_STRING}\n"
    # Use string in sed -i
    echo -ne " * Processing CosmicFusionsList.tsv with sed\n"
    sed -i -e "${SED_STRING}" CosmicFusionsList.tsv
fi

# Remove existing jobs file if exists
if [[ -e Arriba_jobs.txt ]]
then
    rm Arriba_jobs.txt
fi

# Generate Arriba jobs
echo -ne "\nGenerating Arriba jobs...\n"
COUNTER=1
while [[ $COUNTER -le $SAMPLES ]];
do
    echo -ne "  - Working on sample $COUNTER of $SAMPLES, ID: "
    SAMPLE_ID=$(awk "NR==${COUNTER}" ${SAMPLE_SHEET} | cut -f 1)
    echo -ne "${SAMPLE_ID}"
    echo -ne "\n"
    echo "${ARRIBA} ${ARRIBA_ARG} -U ${MAX_READS} -x ${OUTPUT_DIR}/${SAMPLE_ID}.bam -o ${OUTPUT_DIR}/${SAMPLE_ID}_fusions.tsv -O ${OUTPUT_DIR}/${SAMPLE_ID}_fusions.discarded.tsv -a ${REF}.fa -g ${ANNO}.gtf -k CosmicFusionsList.tsv -b ${BLACKLIST} -T -T -P -P" >> Arriba_jobs.txt
    ((COUNTER++))
done

# Generate SAMtools jobs and remove unsorted BAM
# We need to sort and index our BAM so we can view alignments in samtools

# Remove existing jobs file if exists
if [[ -e SAMtools_jobs.txt ]]
then
    rm SAMtools_jobs.txt
fi

echo -ne "\nGenerating SAMtools jobs to sort output BAM...\n"
COUNTER=1
while [[ $COUNTER -le $SAMPLES ]];
do
    echo -ne "  - Working on sample $COUNTER of $SAMPLES, ID: "
    SAMPLE_ID=$(awk "NR==${COUNTER}" ${SAMPLE_SHEET} | cut -f 1)
    echo -ne "${SAMPLE_ID}.bam will be sorted as ${SAMPLE_ID}_Sorted.bam and indexed, gzip compression level is set to ${SAMTOOLS_CMPLVL}\n"
    echo "samtools sort -@ ${SAMTOOLS_CPU} -l ${SAMTOOLS_CMPLVL} -m ${SAMTOOLS_RAM} ${OUTPUT_DIR}/${SAMPLE_ID}.bam  > ${OUTPUT_DIR}/${SAMPLE_ID}_sorted.bam; samtools index ${OUTPUT_DIR}/${SAMPLE_ID}_sorted.bam; rm ${OUTPUT_DIR}/${SAMPLE_ID}.bam" >> SAMtools_jobs.txt
    ((COUNTER++))
done

# Generate plotting jobs in R

# Remove existing jobs file if exists
if [[ -e Plotting_jobs.txt ]]
then
    rm Plotting_jobs.txt
fi

echo -ne "\nGenerating Plotting jobs...\n"
COUNTER=1
while [[ $COUNTER -le $SAMPLES ]];
do
    echo -ne "  - Working on sample $COUNTER of $SAMPLES, ID: "
    SAMPLE_ID=$(awk "NR==${COUNTER}" ${SAMPLE_SHEET} | cut -f 1)
    echo -ne "${SAMPLE_ID}, plots will be written to ${OUTPUT_DIR}/${SAMPLE_ID}.pdf\n"
    echo "${ARRIBA_DIR}/draw_fusions.R --fusions=${OUTPUT_DIR}/${SAMPLE_ID}_fusions.tsv --alignments=${OUTPUT_DIR}/${SAMPLE_ID}_sorted.bam --output=${OUTPUT_DIR}/${SAMPLE_ID}.pdf --annotation=${ANNO}.gtf --cytobands=${CYTOBANDS} --proteinDomains=${DOMAINS}" >> Plotting_jobs.txt
    ((COUNTER++))
done


# Run STAR Jobs
tput bold
echo -ne "\nRunning STAR jobs...\n"
tput sgr0
parallel --progress --jobs ${GNUP_THREAD_STAR} --joblog STAR_joblog.txt < STAR_jobs.txt

# Run Arriba jobs
tput bold
echo -ne "\nRunning Arriba jobs...\n"
tput sgr0
parallel --progress --jobs ${GNUP_THREAD} --joblog Arriba_joblog.txt < Arriba_jobs.txt

# Run SAMtools jobs
tput bold
echo -ne "\nRunning SAMtools jobs...\n"
tput sgr0
parallel --progress --jobs ${GNUP_THREAD_SAMTOOLS} --joblog SAMtools_joblog.txt < SAMtools_jobs.txt

# Run Plotting jobs
tput bold
echo -ne "\nRunning Plotting jobs...\n"
tput sgr0
parallel --progress --jobs ${GNUP_THREAD} --joblog Plotting_joblog.txt < Plotting_jobs.txt

tput bold
echo -ne "\n\nDone!\n\n"
tput sgr0
