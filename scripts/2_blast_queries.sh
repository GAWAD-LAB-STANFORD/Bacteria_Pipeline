#!/bin/bash
#
#SBATCH --job-name=2_blast_queries
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
SCRIPT_COMMAND="$@"
while [ "$1" != "" ]; do
    case $1 in
        --scratch_dir )             shift
                                    SCRATCH_DIR=$1
                                    ;;
        --tools_dir )               shift
                                    TOOLS_DIR=$1
                                    ;;
        --blast_db_types_string )   shift
                                    BLAST_DB_TYPES_ARRAY=( $(echo $1 | sed 's/-/ /g') )
                                    ;;
        --ncbi_db_dir_prefix )      shift
                                    NCBI_DB_DIR_PREFIX=$1
                                    ;;
        --identify )                shift
                                    IDENTIFY=$1
                                    ;;
        --blast_num_alignments )    shift
                                    BLAST_NUM_ALIGNMENTS=$1
                                    ;;
        --blast_align_min )         shift
                                    BLAST_ALIGN_MIN=$1
                                    ;;
        --query )                   shift
                                    QUERY=$1
                                    ;;
        --long_query_string )       shift
                                    LONG_QUERY_ARRAY=( $(echo $1 | sed 's/:/ /g') )
                                    ;;
    esac
    shift
done

if [ -z $SCRATCH_DIR ] || [ -z $TOOLS_DIR ] || [ -z $BLAST_DB_TYPES_ARRAY ] || [ -z $NCBI_DB_DIR_PREFIX ] || \
    [ -z $IDENTIFY ] || [ -z $BLAST_NUM_ALIGNMENTS ] || [ -z $BLAST_ALIGN_MIN ] || [ -z $QUERY ] || \
    [ -z $LONG_QUERY_ARRAY ]; then
    echo "Variables not supplied correctly. Check script for intake parameters. All are required to be specified. Exiting with code 1"
    exit 1
fi

BLAST_QUERY=${LONG_QUERY_ARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}
QUERY_NUM=$(echo $BLAST_QUERY | sed "s/${IDENTIFY}_long_${QUERY}s_//" | sed "s/.fasta//")
echo -e "START: $(date)\nBacteria Pipeline\nScript command: $SCRIPT_COMMAND\n${QUERY}: $BLAST_QUERY"
cd $SCRATCH_DIR

for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
    echo "DB type: $DB_TYPE"
    export BLASTDB=${NCBI_DB_DIR_PREFIX}${DB_TYPE}
    echo "${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $DB_TYPE -num_alignments $BLAST_NUM_ALIGNMENTS -perc_identity $BLAST_ALIGN_MIN \
        -num_threads 4 -outfmt 15 -query $BLAST_QUERY -out ${IDENTIFY}_blast_results_${DB_TYPE}_${QUERY_NUM}.json"
    ${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $DB_TYPE -num_alignments $BLAST_NUM_ALIGNMENTS -perc_identity $BLAST_ALIGN_MIN \
        -num_threads 2 -outfmt 15 -query $BLAST_QUERY -out ${IDENTIFY}_blast_results_${DB_TYPE}_${QUERY_NUM}.json
done
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"