#!/bin/bash
#
#SBATCH --job-name=2_blast_contigs
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
RESULTS_DIR=$1
TOOLS_DIR=$2
BLAST_DB_TYPES_ARRAY=( $(echo $3 | sed 's/-/ /g') )
NCBI_DB_DIR_PREFIX=$4
IDENTIFY=$5
LONG_CONTIG_ARRAY=( $(echo $6 | sed 's/:/ /g') )
CONTIG=${LONG_CONTIG_ARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}
CONTIG_NUM=$(echo $CONTIG | sed "s/${IDENTIFY}_long_contigs_//" | sed "s/.fasta//")

echo -e "START: $(date)\nBacteria Pipeline\nResults dir: $RESULTS_DIR\nSlurm ID: $SLURM_ARRAY_TASK_ID\nContig: $CONTIG"
cd $RESULTS_DIR

for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
    echo "DB type: $DB_TYPE"
    export BLASTDB=${NCBI_DB_DIR_PREFIX}${DB_TYPE}
    ${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $DB_TYPE -num_alignments 5 -num_threads 4 -outfmt 15 -query $CONTIG \
        -out ${IDENTIFY}_blast_results_${DB_TYPE}_${CONTIG_NUM}.json
done
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"