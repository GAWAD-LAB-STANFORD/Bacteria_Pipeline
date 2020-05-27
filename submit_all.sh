#!/bin/bash

PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
HELP="\
Purpose: \n\t\
    This pipeline is built to identify bacterial species from pair-end fastq.gz files and remove human contamination \n\n\
Required arguments: -p/--project <arg> and either -f/--fastq_dir <arg> or -r/--results_dir <arg> \n\
Optional arguments: -b/--run_dir <arg>, --sample_sheet <arg>, --R1_suffix <arg>, --R2_suffix <arg>, --err_out_dir <arg>, --contig_len_min <arg>, --slurm <arg> \n\
Defaults: \n\t\
    If no fastq_dir specified, uses results_dir \n\t\
    If no results_dir specified, makes new directory in fastq_dir \n\t\
    sample_sheet: SampleSheet.csv \n\t\
    R1_suffix: _L001_R1_001.fastq.gz or _R1_001.fastq.gz \n\t\
    R2_suffix: _L001_R2_001.fastq.gz or _R1_001.fastq.gz \n\t\
    min_contig_len: 5000 \n\n\
Run after demultiplexing: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --fastq_dir /oak/stanford/groups/cgawad/MRD_project/ --project MRD_project \n\n\
Run with demultiplexing: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/191126_MN01236_0003_A000H2WWHT --fastq_dir /oak/stanford/groups/cgawad/MRD_project/ --project MRD_project \n\n\
Run with demultiplexing, wait 12 hours before starting, and email notification when analysis begins and ends: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/191126_MN01236_0003_A000H2WWHT --fastq_dir /oak/stanford/groups/cgawad/MRD_project/ --project MRD_project --slurm --begin=now+12hours --mail-type=ALL \n\n\
For more information, read the README.md"

# Reads in command line option arguments and assigns them to variables
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )           echo -e $HELP
                                exit 0
                                ;;
        -b | --run_dir )        shift
                                RUN_DIR=$1
                                ;;
        --sample_sheet )        shift
                                SAMPLE_SHEET=$1
                                ;;
        --R1_suffix )           shift
                                R1_SUFFIX=$1
                                ;;
        --R2_suffix )           shift
                                R2_SUFFIX=$1
                                ;;
        --err_out_dir )         shift
                                STD_ERR_OUT_DIR=$1
                                ;;
        --contig_len_min )      shift
                                CONTIG_LENGTH_MINIMUM=$1
                                ;;
        -f | --fastq_dir )      shift
                                FASTQ_DIR=$1
                                ;;
        -r | --results_dir )    shift
                                RESULTS_DIR=$1
                                ;;
        -p | --project )        shift
                                PROJECT=$1
                                ;;
        --slurm )               shift
                                SLURM_OPTIONS=${@:1}
                                ;;
    esac
    shift
done

# Ensure we have the requires variables set and set other variables
if ([ -z $FASTQ_DIR ] && [ -z $RESULTS_DIR ]) || [ -z $PROJECT ] || [ -z $PIPELINE_DIR ]; then
    echo "Variables not supplied correctly. Use -h/--help options for assistance. Ending program..."
    exit 1
fi
if [ -z $FASTQ_DIR ]; then
    FASTQ_DIR="$RESULTS_DIR"
elif [ -z $RESULTS_DIR ]; then
    RESULTS_DIR="${FASTQ_DIR}/$(date '+%Y-%m-%d')_${PROJECT}_Results"
fi
if [ -z $STD_ERR_OUT_DIR ]; then
    STD_ERR_OUT_DIR="${RESULTS_DIR}/std_err_out_files"
fi
# Make results and std error output directories if they don't exist
if [ ! -d $FASTQ_DIR ]; then
    mkdir $FASTQ_DIR
fi
if [ ! -d $RESULTS_DIR ]; then
    mkdir $RESULTS_DIR
fi
if [ ! -d $STD_ERR_OUT_DIR ]; then
    mkdir $STD_ERR_OUT_DIR
fi
OPTIONS=()
if [ ! -z $RUN_DIR ] && [ ! -z $SAMPLE_SHEET ]; then
    OPTIONS+=( "--run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET" )
elif [ ! -z $RUN_DIR ]; then
    OPTIONS+=( "--run_dir $RUN_DIR" )
elif [ ! -z $SAMPLE_SHEET ]; then
    echo "Variables not supplied correctly. Please specify a run diretory for demultiplexing with --run_dir. Exiting with code 1"
    exit 1
fi
if [ ! -z $R1_SUFFIX ]; then
    OPTIONS+=( "--R1_suffix $R1_SUFFIX" )
fi
if [ ! -z $R2_SUFFIX ]; then
    OPTIONS+=( "--R2_suffix $R2_SUFFIX" )
fi
if [ ! -z $CONTIG_LENGTH_MINIMUM ]; then
    OPTIONS+=( "--contig_len_min $CONTIG_LENGTH_MINIMUM" )
fi

sbatch ${SLURM_OPTIONS[@]} -J $PROJECT -e $STD_ERR_OUT_DIR/%A_%x.err -o $STD_ERR_OUT_DIR/%A_%x.out ${PIPELINE_DIR}/pipeline_control.sh --err_out_dir $STD_ERR_OUT_DIR -f $FASTQ_DIR -r $RESULTS_DIR -p $PROJECT -d $PIPELINE_DIR ${OPTIONS[@]}