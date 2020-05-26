#!/bin/bash
#
#SBATCH --job-name=demultiplexer
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=2
#SBATCH --time=6:00:00
#SBATCH --partition=cgawad

while [ "$1" != "" ]; do
    case $1 in
        --run_dir )         shift
                            RUN_DIR=$1
                            ;;
        --sample_sheet )    shift
                            SAMPLE_SHEET=$1
                            ;;
        --fastq_dir )       shift
                            FASTQ_DIR=$1
                            ;;
    esac
    shift
done

if [ -z $RUN_DIR ] || [ -z $FASTQ_DIR ]; then
    echo "Variables not supplied correctly. Use -h/--help options for assistance. Exiting with code 1"
    exit 1
fi
if [ -z $SAMPLE_SHEET ]; then
    SAMPLE_SHEET="${RUN_DIR}/SampleSheet.csv"
fi
if [ ! -f $SAMPLE_SHEET ]; then
    echo "Sample sheet $SAMPLE_SHEET not found. Exiting with code 1"
    exit 1
fi

ml biology bcl2fastq
bcl2fastq --runfolder-dir $RUN_DIR --sample-sheet $SAMPLE_SHEET --output-dir $FASTQ_DIR
