#!/bin/bash
#
#SBATCH --job-name=2_blast_queries
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
SCRATCH_DIR=$1
TOOLS_DIR=$2
BLAST_DB_TYPES_ARRAY=( $(echo $3 | sed 's/-/ /g') )
NCBI_DB_DIR_PREFIX=$4
IDENTIFY=$5
BLAST_NUM_ALIGNMENTS=$6
BLAST_ALIGN_MIN=$7
QUERY=$8
LONG_QUERY_ARRAY=( $(echo $9 | sed 's/:/ /g') )
BLAST_QUERY=${LONG_QUERY_ARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}
QUERY_NUM=$(echo $BLAST_QUERY | sed "s/${IDENTIFY}_long_${QUERY}s_//" | sed "s/.fasta//")

echo -e "START: $(date)\nBacteria Pipeline\nScratch dir: $SCRATCH_DIR\n${QUERY}: $BLAST_QUERY"
cd $SCRATCH_DIR

for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
    echo "DB type: $DB_TYPE"
    export BLASTDB=${NCBI_DB_DIR_PREFIX}${DB_TYPE}
    echo "${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $DB_TYPE -num_alignments $BLAST_NUM_ALIGNMENTS -perc_identity $BLAST_ALIGN_MIN \
        -num_threads 4 -outfmt 15 -query $BLAST_QUERY -out ${IDENTIFY}_blast_results_${DB_TYPE}_${QUERY_NUM}.json"
    ${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $DB_TYPE -num_alignments $BLAST_NUM_ALIGNMENTS -perc_identity $BLAST_ALIGN_MIN \
        -num_threads 4 -outfmt 15 -query $BLAST_QUERY -out ${IDENTIFY}_blast_results_${DB_TYPE}_${QUERY_NUM}.json
done
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"