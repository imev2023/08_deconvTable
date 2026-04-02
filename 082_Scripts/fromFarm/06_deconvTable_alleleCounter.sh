#!/bin/bash

#BSUB -M 10GB # Memory
#BSUB -R "select[mem>10GB] rusage[mem=10GB] span[hosts=1]"
#BSUB -n 8 # CPU's - matches hosts
#BSUB -q normal # Queue (n=normal, l=long)
#BSUB -J "alleleCounter[1-2]"
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/063_Results/06_deconvTable_logs/06_deconvTable_alleleCounter_qNormal_10GB.%J.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/063_Results/06_deconvTable_logs/06_deconvTable_alleleCounter_qNormal_10GB.%J.err

# run this from the output dir bsub -G team267-grp < /lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/062_Scripts/06_deconvTable_alleleCounter.sh

# load modules
module load ISG/singularity/3.11.4
module load alleleCount/4.3.0

# set environmental variables
WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/"
DATADIR="$WDIR/061_Data/"
CRAMDIR="$DATADIR/061_crams/"
RESDIR="$WDIR/063_Results/06_deconvTable_alleleCounter/"
REF_IDX="/lustre/scratch125/casm/teams/team267_murchison/ref_genomes/Sarcophilus_Harrisii/Ensembl_mSarHar1.11/Sarcophilus_harrisii.mSarHar1.11.dna.toplevel.fa.fai"
LOGDIR="$WDIR/063_Results/06_deconvTable_logs/"
mkdir -p $RESDIR

# Define the array of CRAM files
FILES=("$CRAMDIR"/*.cram)            # List all CRAM files in the directory
NUM_FILES=${#FILES[@]}                # Count the number of files

# Get the current file based on LSB_JOBINDEX
FILE_INDEX=$((LSB_JOBINDEX - 1))           # Adjust to zero-indexed array
FILE="${FILES[$FILE_INDEX]}"          # Select the corresponding file

# Define output file name
OUTPUT=$(basename "${FILE%.cram}")         

# Create log folder per sample
LOG_DIR="$LOGDIR/$OUTPUT" 
mkdir -p "$LOG_DIR"

# create ouput folder and cd
alleleCounter \
    -l "$DATADIR/061_positions/08_deconvTable_dft2SomRefAltExport_withoutColnames.tsv" \
    -b "$FILE" \
    -o "$RESDIR/06_deconvTable_dft2SomRefAltExport_alleleCounter_${OUTPUT}" \
    --ref-file "$REF_IDX" 