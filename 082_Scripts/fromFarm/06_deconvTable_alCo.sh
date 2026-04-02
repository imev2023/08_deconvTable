#!/bin/bash

#BSUB -M 5GB # Memory
#BSUB -R "select[mem>5GB] rusage[mem=5GB] span[hosts=1]"
#BSUB -n 1 # CPU's - matches hosts
#BSUB -q normal # Queue (n=normal, l=long)
#BSUB -J "dft1Som[66]"
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/063_Results/06_deconvTable_logs/06_deconvTable_alCo_qNormal_dft1Som_5GB.%J.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/063_Results/06_deconvTable_logs/06_deconvTable_alCo_qNormal_dft1Som_5GB.%J.err

# run this from the output dir bsub -G team267-grp < /lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/062_Scripts/06_deconvTable_alCo.sh

# load modules
module load ISG/singularity/3.11.4

# set environmental variables
WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/06_deconvTable/"
DATADIR="$WDIR/061_Data/"
CRAMDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/05_devilRna/053_Results/053_devilRna_crams/0531_devilRna_maxRna_realigned/"
RESDIR="$WDIR/063_Results/06_deconvTable_alCo_dft1SomRefAltExport/"
REF_IDX="/lustre/scratch125/casm/teams/team267_murchison/ref_genomes/Sarcophilus_Harrisii/Ensembl_mSarHar1.11/Sarcophilus_harrisii.mSarHar1.11.dna.toplevel.fa.fai"
LOGDIR="$WDIR/063_Results/06_deconvTable_logs/"
mkdir -p $RESDIR

# Define the array of CRAM files
FILES=("$CRAMDIR"/*T*.cram)            # List all CRAM files in the directory
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
echo "Processing $FILE : $(date)" > "$LOG_DIR/06_deconvTable_dft1SomRefAltExport_alCo_${OUTPUT}_qNormal_5GB.log"

singularity run -B '/nfs:/nfs' -B '/lustre:/lustre' /nfs/dog_n_devil/kevin/code/rust/alco/alco.sif \
    --bamfile "$FILE" \
    --locifile "$DATADIR/061_positions/08_deconvTable_dft1SomRefAltExport_withoutColnames.tsv" \
    > "$RESDIR/06_deconvTable_dft1SomRefAltExport_${OUTPUT}.txt" \
    2> "$LOG_DIR/06_deconvTable_${OUTPUT}.err"

echo "End alCo $FILE : $(date)" >> "$LOG_DIR/06_deconvTable_dft1SomRefAltExport_alCo_${OUTPUT}_qNormal_5GB.log"