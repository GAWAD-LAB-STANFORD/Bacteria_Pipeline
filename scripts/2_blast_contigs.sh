#!/bin/bash
#
#SBATCH --job-name=2_blast_contigs
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
RESULTS_DIR=$1
ORGANIZED_CONTIGS_PREFIX=$2
NCBI_BLAST_TOOL_DIR=$3
BLAST_DB_TYPES_ARRAY=( $(echo $4 | sed 's/-/ /g') )
NCBI_DB_DIR_PREFIX=$5
BLAST_RESULTS_PREFIX=$6

echo -e "START: $(date)\nBacteria Pipeline\nResults dir: $RESULTS_DIR\nSlurm ID: $SLURM_ARRAY_TASK_ID"
cd $RESULTS_DIR
CONTIG=$(ls ${ORGANIZED_CONTIGS_PREFIX}* | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Contig: $CONTIG"
for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
    echo "DB type: $DB_TYPE"
    export BLASTDB=${NCBI_DB_DIR_PREFIX}${DB_TYPE}
    ${NCBI_BLAST_TOOL_DIR}/bin/blastn -db $DB_TYPE -num_alignments 5 -num_threads 4 -outfmt 15 -query $CONTIG \
        -out ${BLAST_RESULTS_PREFIX}${DB_TYPE}_${SLURM_ARRAY_TASK_ID}.json
done
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"