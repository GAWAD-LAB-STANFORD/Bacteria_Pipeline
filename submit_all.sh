#!/bin/bash
#
#SBATCH --job-name=submit_all
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=5:00:00
#SBATCH --partition=cgawad

PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
PIPELINE_COMMAND="$@"
HELP="\
Purpose: \n\t\
    To identify bacterial species from pair-end fastq.gz files and remove human contamination \n\n\
Required arguments: -p/--project <arg> and either -f/--fastq_dir <arg> or -r/--results_dir <arg> \n\
Optional arguments: -s/--scratch_dir <arg>, --err_out_dir <arg>, --skip_scratch, -b/--run_dir <arg>, \n\t\
    --sample_sheet <arg>, --skip_identify, --only_identify, --identify <arg>, \n\t\
    --R1_suffix <arg>, --R2_suffix <arg>, --skip_trimming, --rna, --filter_rhesus, \n\t\
    --contig_len_min <arg>, --kraken_db_types <arg>, --blast_db_types <arg>, \n\t\
    --num_alignments <arg>, --contig_align_min <arg>, --add_genus, --slurm <arg> \n\
Defaults: \n\t\
    If no fastq_dir specified, uses results_dir \n\t\
    If no results_dir specified, makes new directory in fastq_dir \n\t\
    scratch_dir: /scratch/groups/cgawad/date_project_Scratch \n\t\
    sample_sheet: SampleSheet.csv \n\t\
    R1_suffix: _L001_R1_001.fastq.gz or _R1_001.fastq.gz or _R1.fastq.gz \n\t\
    R2_suffix: _L001_R2_001.fastq.gz or _R2_001.fastq.gz or _R2.fastq.gz \n\t\
    contig_len_min: 5000 \n\t\
    kraken_db_types: microbial \n\t\
    blast_db_types: nt \n\t\
    num_alignments: 5 \n\t\
    contig_align_min: 0.9 \n\n\
Run after demultiplexing and with fastq directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --project 2020-01-01_Project \n\n\
Run after demultiplexing and with results directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --results_dir /oak/stanford/groups/cgawad/2020-01-01_Results/ --project 2020-01-01_Project \n\n\
Run with demultiplexing and fastq directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/2020-01-01_BCLs --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --project 2020-01-01_Project \n\n\
Run with demultiplexing and results directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/2020-01-01_BCLs --results_dir /oak/stanford/groups/cgawad/2020-01-01_Results/ --project 2020-01-01_Project \n\n\
Run with demultiplexing, fastq directory, and results directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/2020-01-01_BCLs --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --results_dir /oak/stanford/groups/cgawad/2020-01-01_Results/ --project 2020-01-01_Project \n\n\
For more information, read the README.md"

# Reads in command line option arguments and assigns them to variables
SKIP_SCRATCH=0
SKIP_IDENTIFY=0
ONLY_IDENTIFY=0
SKIP_TRIMMOMATIC=0
RNA=0
FILTER_RHESUS=0
CONTIG_LENGTH_MINIMUM=5000
KRAKEN_DB_TYPES="microbial"
BLAST_DB_TYPES="nt"
NUM_ALIGNMENTS=5
CONTIG_ALIGN_MINIMUM=0.9
ADD_GENUS=0
STEP=0
TEMP_ARRAY_START=0
DEPENDENCY=""
DEPENDER=""
FIGURE_OPTIONS=()
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )           echo -e $HELP
                                exit 0
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
        -d | --pipeline_dir )   shift
                                PIPELINE_DIR=$1
                                ;;
        -s | --scratch_dir )    shift
                                SCRATCH_DIR=$1
                                ;;
        --err_out_dir )         shift
                                STD_ERR_OUT_DIR=$1
                                ;;
        --skip_scratch )        SKIP_SCRATCH=1
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
        --skip_identify )       SKIP_IDENTIFY=1
                                ;;
        --only_identify )       ONLY_IDENTIFY=1
                                ;;
        --identify )            shift
                                IDENTIFY=$1
                                ;;
        --skip_trimming )       SKIP_TRIMMOMATIC=1
                                ;;
        --rna )                 RNA=1
                                ;;
        --filter_rhesus )       FILTER_RHESUS=1
                                ;;
        --contig_len_min )      shift
                                CONTIG_LENGTH_MINIMUM=$1
                                ;;
        --add_genus )           shift
                                ADD_GENUS=1
                                ;;
        --kraken_db_types )     shift
                                KRAKEN_DB_TYPES=$1
                                ;;
        --blast_db_types )      shift
                                BLAST_DB_TYPES=$1
                                ;;
        --num_alignments )      shift
                                NUM_ALIGNMENTS=$1
                                ;;
        --contig_align_min )    shift
                                CONTIG_ALIGN_MIN=$1
                                ;;
        --step1 )               STEP=1
                                ;;
        --step2 )               STEP=2
                                ;;
        --step3 )               STEP=3
                                ;;
        --slurm )               shift
                                SLURM_OPTIONS=${@:1}
                                ;;
    esac
    shift
done

# Hardcoded paths and variables
TEMP_ARRAY_INCREMENT=1000
REFERENCE_DIR="/oak/stanford/groups/cgawad/Reference_Files/GATK_Resource_Bundle_hg38"
TOOLS_DIR="/oak/stanford/groups/cgawad/Sequencing_Analysis_Tools"
REF_FASTA="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
RHESUS_FASTA="/oak/stanford/groups/cgawad/Reference_Files/Macaca_mulattta_Rhesus_monkey_hg38/Macaca_mulatta_Rhesus_monkey_hg38.fasta"
SCRIPT_DIR="${PIPELINE_DIR}/scripts"
KRAKEN_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/Kraken2_Fatfree_Databases/kraken2-fatfree-"
NCBI_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/NCBI_RefSeq_Databases/ncbi_database_"
NCBI_ANNOTATIONS_DIR="/oak/stanford/groups/cgawad/Reference_Files/NCBI_Annotations"
CONTIGS_PER_BLAST_JOB=320

# Ensure we have the required variables set and set other variables
if ([ -z $FASTQ_DIR ] && [ -z $RESULTS_DIR ]) || [ -z $PROJECT ] || [ -z $PIPELINE_DIR ]; then
    echo "Variables not supplied correctly. Use -h/--help options for assistance. Ending program..."
    exit 1
fi
if [ -z $FASTQ_DIR ]; then
    FASTQ_DIR="$RESULTS_DIR"
elif [ -z $RESULTS_DIR ]; then
    RESULTS_DIR="${FASTQ_DIR}/$(date '+%Y-%m-%d')_${PROJECT}_Results"
fi
if [ -z $SCRATCH_DIR ] && [ $SKIP_SCRATCH -eq 0 ]; then
    SCRATCH_DIR="/scratch/groups/cgawad/$(date '+%Y-%m-%d')_${PROJECT}_Scratch"
elif [ $SKIP_SCRATCH -eq 1 ]; then
    SCRATCH_DIR="$RESULTS_DIR"
fi
if [ -z $STD_ERR_OUT_DIR ]; then
    STD_ERR_OUT_DIR="${RESULTS_DIR}/std_err_out_files"
fi
# Make directories if they don't exist
if [ ! -d $FASTQ_DIR ]; then
    mkdir $FASTQ_DIR
fi
if [ ! -d $RESULTS_DIR ]; then
    mkdir $RESULTS_DIR
fi
if [ ! -d $SCRATCH_DIR ]; then
    mkdir $SCRATCH_DIR
fi
if [ ! -d $STD_ERR_OUT_DIR ]; then
    mkdir $STD_ERR_OUT_DIR
fi
OPTIONS=( "-f $FASTQ_DIR -r $RESULTS_DIR -d $PIPELINE_DIR -p $PROJECT -s $SCRATCH_DIR --err_out_dir $STD_ERR_OUT_DIR " )
if [ ! -z $RUN_DIR ] && [ $ONLY_IDENTIFY -eq 1 ]; then
    echo "Variables not supplied correctly. Cannot perform demultiplexing while only identifying data. Exiting with code 1"
    exit 1
fi
if [ ! -z $RUN_DIR ] && [ -z $SAMPLE_SHEET ]; then
    SAMPLE_SHEET="${RUN_DIR}/SampleSheet.csv"
elif [ -z $RUN_DIR ] && [ ! -z $SAMPLE_SHEET ]; then
    echo "Variables not supplied correctly. Please specify a run diretory for demultiplexing with --run_dir. Exiting with code 1"
    exit 1
fi
if [ ! -z $RUN_DIR ] && [ ! -z $SAMPLE_SHEET ]; then
    if [ ! -f $SAMPLE_SHEET ]; then
        echo "Sample sheet $SAMPLE_SHEET not found. Exiting with code 1"
        exit 1
    fi
    OPTIONS+=( "--run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET" )
fi
if [ $SKIP_IDENTIFY -eq 1 ] && [ $ONLY_IDENTIFY -eq 1 ]; then
    echo "Variables not supplied correctly. Please specify either --skip_identify or --only_identify, not both. Exiting with code 1"
    exit 1
elif [ $SKIP_IDENTIFY -eq 1 ]; then
    OPTIONS+=( "--skip_identify" )
elif [ $ONLY_IDENTIFY -eq 1 ]; then
    OPTIONS+=( "--only_identify" )
fi
if [ ! -z $IDENTIFY ]; then
    OPTIONS+=( "--identify $IDENTIFY" )
else
    IDENTIFY=$PROJECT
fi
if [ ! -z $R1_SUFFIX ]; then
    OPTIONS+=( "--R1_suffix $R1_SUFFIX" )
fi
if [ ! -z $R2_SUFFIX ]; then
    OPTIONS+=( "--R2_suffix $R2_SUFFIX" )
fi
if [ $SKIP_TRIMMOMATIC -eq 1 ]; then
    OPTIONS+=( "--skip_trimming" )
    SCRATCH_DIR=$RESULTS_DIR
fi
if [ $RNA -eq 1 ]; then
    OPTIONS+=( "--rna" )
fi
if [ $FILTER_RHESUS -eq 1 ]; then
    OPTIONS+=( "--filter_rhesus" )
    REF_FASTA_STRING="${REF_FASTA}:${RHESUS_FASTA}"
    REF_NAME_STRING="human:rhesus"
else
    REF_FASTA_STRING="${REF_FASTA}"
    REF_NAME_STRING="human"
fi
if [ $CONTIG_LENGTH_MINIMUM -ne 5000 ]; then
    OPTIONS+=( "--contig_len_min $CONTIG_LENGTH_MINIMUM" )
fi
if [ "$KRAKEN_DB_TYPES" != "microbial" ]; then
    OPTIONS+=( "--kraken_db_types $KRAKEN_DB_TYPES" )
fi
if [ "$BLAST_DB_TYPES" != "nt" ]; then
    OPTIONS+=( "--blast_db_types $BLAST_DB_TYPES" )
fi
if [ $NUM_ALIGNMENTS -ne 250 ]; then
    OPTIONS+=( "--num_alignments $NUM_ALIGNMENTS" )
fi
if [ "$CONTIG_ALIGN_MINIMUM" = "0.9" ]; then
    OPTIONS+=( "--contig_align_min $CONTIG_ALIGN_MINIMUM" )
fi
if [ $ADD_GENUS -eq 1 ]; then
    FIGURE_OPTIONS+=( "--add_genus" )
fi


TEMP_PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
if [ $ONLY_IDENTIFY -eq 1 ]; then
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${IDENTIFY}_pipeline_status.txt
else
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${PROJECT}_pipeline_status.txt
fi
cd $SCRATCH_DIR
if [ "$TEMP_PIPELINE_DIR" = "$PIPELINE_DIR" ]; then
    echo -e "\nSTART: $(date)\nBacteria Pipeline\n\n$PIPELINE_COMMAND\n\nProject: $PROJECT\nResults dir: $RESULTS_DIR\nFastq dir: $FASTQ_DIR\nScratch dir: $SCRATCH_DIR\nErr out dir: $STD_ERR_OUT_DIR" >> $PIPELINE_STATUS
    # Optional variable definition
    if [ $SKIP_SCRATCH -eq 0 ]; then
        echo "Default: Scratch dir is different from Results dir" >> $PIPELINE_STATUS
    else
        echo "Option: Scratch dir is the same as Results dir" >> $PIPELINE_STATUS
    fi
    if [ $SKIP_IDENTIFY -eq 1 ]; then
        echo "Option: Skip identification of data - will only process the fastqs, build the contigs, and run Kraken2" >> $PIPELINE_STATUS
    fi
    if [ $ONLY_IDENTIFY -eq 1 ]; then
        echo "Option: Only identification of data - will only BLAST and filter from already built contigs" >> $PIPELINE_STATUS
    fi
    if [ "$IDENTIFY" != "$PROJECT" ]; then
        echo "Option: Identify different from Project: $IDENTIFY" >> $PIPELINE_STATUS
    fi
    if [ $SKIP_TRIMMOMATIC -eq 1 ]; then
        echo "Option: Skip trimming - will not run trimmomatic" >> $PIPELINE_STATUS
    fi
    if [ $RNA -eq 1 ]; then
        echo "Option: Expecting RNA input and will align using STAR instead of BWA" >> $PIPELINE_STATUS
    fi
    if [ $FILTER_RHESUS -eq 1 ]; then
        echo "Option: Filter rhesus - will remove reads that align to macaca mulatta rhesus monkey" >> $PIPELINE_STATUS
    fi
    if [ $CONTIG_LENGTH_MINIMUM -eq 5000 ]; then
        echo "Default: Contig length minimum: 5000" >> $PIPELINE_STATUS
    else
        echo "Option: Contig length minimum: $CONTIG_LENGTH_MINIMUM" >> $PIPELINE_STATUS
    fi
    if [ "$KRAKEN_DB_TYPES" = "microbial" ]; then
        echo "Default: Kraken db types: microbial" >> $PIPELINE_STATUS
    else
        echo "Option: Kraken db types: $KRAKEN_DB_TYPES" >> $PIPELINE_STATUS
    fi
    if [ "$BLAST_DB_TYPES" = "nt" ]; then
        echo "Default: Blast db types: nt" >> $PIPELINE_STATUS
    else
        echo "Option: Blast db types: $BLAST_DB_TYPES" >> $PIPELINE_STATUS
    fi
    if [ "$CONTIG_ALIGN_MINIMUM" = "0.9" ]; then
        echo "Default: Contig align minimum: 0.9" >> $PIPELINE_STATUS
    else
        echo "Option: Contig align minimum: $CONTIG_ALIGN_MINIMUM" >> $PIPELINE_STATUS
    fi
    if [ $ADD_GENUS -eq 1 ]; then
        echo "Option: Adding genus to figures" >> $PIPELINE_STATUS
    fi
    echo " " >> $PIPELINE_STATUS
fi


if [ ! -z $SLURM_OPTIONS ]; then
    echo "Option: Slurm - entire pipeline run will be queued with user parameters" >> $PIPELINE_STATUS
    sbatch -J $PROJECT ${SLURM_OPTIONS[@]} \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh ${OPTIONS[@]}
    exit 0
fi


if ([ $STEP -eq 0 ] && [ -z $RUN_DIR ]) || [ $STEP -eq 1 ] || [ $STEP -eq 2 ]; then
    if [ -z $R1_SUFFIX ] || [ -z $R2_SUFFIX ]; then
        R1_SUFFIX="_L001_R1_001.fastq.gz"
        R2_SUFFIX="_L001_R2_001.fastq.gz"
        if [ $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" | wc -l) -eq 0 ]; then
            R1_SUFFIX="_R1_001.fastq.gz"
            R2_SUFFIX="_R2_001.fastq.gz"
        fi
        if [ $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" | wc -l) -eq 0 ]; then
            R1_SUFFIX="_R1.fastq.gz"
            R2_SUFFIX="_R2.fastq.gz"
        fi
    fi
    SAMPLE_ARRAY=( $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" -exec basename {} \; | \
        grep -v "Undetermined" | sed "s/${R1_SUFFIX}//") )
    if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
        echo "No fastq.gz files found in the fastq directory. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    fi
fi


if [ $STEP -eq 0 ] && [ ! -z $RUN_DIR ] && [ $ONLY_IDENTIFY -eq 0 ]; then
    echo "### Demultiplexing ### - START: $(date)" >> $PIPELINE_STATUS
    echo -e "Run dir: $RUN_DIR\nSample sheet: $SAMPLE_SHEET" >> $PIPELINE_STATUS
    echo -e "\nsbatch --parsable -e ${STD_ERR_OUT_DIR}/%A_%x.err -o ${STD_ERR_OUT_DIR}/%A_%x.out \
        ${SCRIPT_DIR}/0_demultiplexer.sh --run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET --fastq_dir $FASTQ_DIR \
        --pipeline_status $PIPELINE_STATUS\n" >> $PIPELINE_STATUS
    DEPENDENCY=$(sbatch --parsable -e ${STD_ERR_OUT_DIR}/%A_%x.err -o ${STD_ERR_OUT_DIR}/%A_%x.out \
        ${SCRIPT_DIR}/0_demultiplexer.sh --run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET --fastq_dir $FASTQ_DIR \
        --pipeline_status $PIPELINE_STATUS)
    echo -e "\nsbatch --dependency=afterok:$DEPENDENCY -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step1 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
    DEPENDER=$(sbatch --dependency=afterok:$DEPENDENCY -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step1 ${OPTIONS[@]})
    echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
elif ([ $STEP -eq 0 ] && [ -z $RUN_DIR ] && [ $ONLY_IDENTIFY -eq 0 ]) || ([ $STEP -eq 1 ] && [ $ONLY_IDENTIFY -eq 0 ]); then
    if [ $TEMP_ARRAY_START -eq 0 ]; then
        echo -e "Number of samples: ${#SAMPLE_ARRAY[@]}\nSamples: ${SAMPLE_ARRAY[@]}\n" >> $PIPELINE_STATUS
        echo "### De novo assembling contigs and detecting contamination ### - START: $(date)" >> $PIPELINE_STATUS
        JOB_COUNT=${#SAMPLE_ARRAY[@]}
        echo "Jobs to run: $JOB_COUNT" >> $PIPELINE_STATUS
        TEMP_ARRAY_START=1
    fi
    
    TEMP_SAMPLE_ARRAY=( ${SAMPLE_ARRAY[@]:$(($TEMP_ARRAY_START - 1)):$TEMP_ARRAY_INCREMENT} )
    TEMP_JOB_COUNT=${#TEMP_SAMPLE_ARRAY[@]}
    echo "Submitting $TEMP_JOB_COUNT jobs for samples $TEMP_ARRAY_START to $(($TEMP_ARRAY_START + ${#TEMP_SAMPLE_ARRAY[@]} - 1))" >> $PIPELINE_STATUS
    TEMP_SAMPLES_STRING=$( IFS=$':'; echo "${TEMP_SAMPLE_ARRAY[*]}" )
    echo -e "\nsbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/1_process_sample.sh \
        $FASTQ_DIR $SCRATCH_DIR $R1_SUFFIX $R2_SUFFIX $SKIP_TRIMMOMATIC $RNA $REF_FASTA_STRING $REF_NAME_STRING \
        $KRAKEN_DB_TYPES $KRAKEN_DB_DIR_PREFIX $TOOLS_DIR $TEMP_SAMPLES_STRING\n" >> $PIPELINE_STATUS
    DEPENDENCY=$(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/1_process_sample.sh \
        $FASTQ_DIR $SCRATCH_DIR $R1_SUFFIX $R2_SUFFIX $SKIP_TRIMMOMATIC $RNA $REF_FASTA_STRING $REF_NAME_STRING \
        $KRAKEN_DB_TYPES $KRAKEN_DB_DIR_PREFIX $TOOLS_DIR $TEMP_SAMPLES_STRING)
    TEMP_ARRAY_START=$(($TEMP_ARRAY_START + $TEMP_ARRAY_INCREMENT))
    echo -e "$(date)\nIncrement: $TEMP_ARRAY_INCREMENT\nNew start: $TEMP_ARRAY_START" >> $PIPELINE_STATUS
    
    if [ $TEMP_ARRAY_START -le ${#SAMPLE_ARRAY[@]} ]; then
        echo -e "\nsbatch --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step1 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        DEPENDER=$(sbatch --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step1 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]})
        echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
    else
        echo -e "\nsbatch --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        DEPENDER=$(sbatch --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 ${OPTIONS[@]})
        echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
    fi
elif [ $STEP -eq 2 ] && [ $ONLY_IDENTIFY -eq 0 ]; then
    SAMPLE_COUNT=1
    for SAMPLE in ${SAMPLE_ARRAY[@]}; do
        if [ ! -f ${SAMPLE}_contigs.fasta ]; then
            echo -e "\tSample number $SAMPLE_COUNT - ${SAMPLE}_contigs.fasta file not found" >> $PIPELINE_STATUS
        else
            rm ${STD_ERR_OUT_DIR}/*_${SAMPLE_COUNT}_1_process_sample.out ${STD_ERR_OUT_DIR}/*_${SAMPLE_COUNT}_1_process_sample.err
        fi
        SAMPLE_COUNT=$((SAMPLE_COUNT+1))
    done
    CONTIG_COUNT=$(ls *_contigs.fasta | wc -l)
    if [ $CONTIG_COUNT -eq 0 ]; then
        echo "No contigs found. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    else
        echo "$CONTIG_COUNT contigs out of a possible ${#SAMPLE_ARRAY[@]} maximum" >> $PIPELINE_STATUS
    fi
    echo "### De novo assembling contigs and detecting contamination ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    ml R/4.2.0
    export R_LIBS="/home/groups/cgawad/R_libs"
    bash ${SCRIPT_DIR}/merge_metrics.sh $PIPELINE_STATUS $PROJECT $KRAKEN_DB_TYPES $RUN_DIR $SAMPLE_SHEET $REF_NAME_STRING
    
    
    if [ $SKIP_IDENTIFY -eq 1 ]; then
        echo "Ending without identfication of data" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 0
    fi
fi


if ([ $STEP -eq 0 ] && [ $ONLY_IDENTIFY -eq 1 ]) || [ $STEP -eq 2 ]; then
    SAMPLE_ARRAY=( $(ls *_contigs.fasta | sed "s/_contigs.fasta//") )
    if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
        echo "No contig fasta files found in the results directory. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    fi
fi


if ([ $STEP -eq 0 ] && [ $ONLY_IDENTIFY -eq 1 ]) || [ $STEP -eq 2 ]; then
    if [ $TEMP_ARRAY_START -eq 0 ]; then
        ml python/3.6.1 biology py-biopython/1.79_py39


        echo "### Organzing contigs ### - START: $(date)" >> $PIPELINE_STATUS
        echo -e "\npython3 ${SCRIPT_DIR}/organize_contigs.py \
            _contigs.fasta $CONTIG_LENGTH_MINIMUM $CONTIGS_PER_BLAST_JOB ${IDENTIFY}_long_contigs_ $PIPELINE_STATUS"
        python3 ${SCRIPT_DIR}/organize_contigs.py \
            "_contigs.fasta" $CONTIG_LENGTH_MINIMUM $CONTIGS_PER_BLAST_JOB "${IDENTIFY}_long_contigs_" $PIPELINE_STATUS
        if [ $(ls ${IDENTIFY}_long_contigs_* | wc -l) -eq 0 ]; then
            echo "No contigs longer than $CONTIG_LENGTH_MINIMUM. Exiting with code 1" >> $PIPELINE_STATUS
            echo "END: $(date)" >> $PIPELINE_STATUS
            exit 1
        fi
        echo "### Organzing contigs ### - END: $(date)" >> $PIPELINE_STATUS
        

        echo "### BLAST aligning contigs ### - START: $(date)" >> $PIPELINE_STATUS
        LONG_CONTIG_ARRAY=( $(ls ${IDENTIFY}_long_contigs_*) )
        JOB_COUNT=${#LONG_CONTIG_ARRAY[@]}
        echo "Jobs to run: $JOB_COUNT" >> $PIPELINE_STATUS
        TEMP_ARRAY_START=1
    fi


    TEMP_LONG_CONTIG_ARRAY=( ${LONG_CONTIG_ARRAY[@]:$(($TEMP_ARRAY_START - 1)):$TEMP_ARRAY_INCREMENT} )
    TEMP_JOB_COUNT=${#TEMP_LONG_CONTIG_ARRAY[@]}
    echo "Submitting $TEMP_JOB_COUNT jobs for samples $TEMP_ARRAY_START to $(($TEMP_ARRAY_START + ${#TEMP_LONG_CONTIG_ARRAY[@]} - 1))" >> $PIPELINE_STATUS
    TEMP_LONG_CONTIGS_STRING=$( IFS=$':'; echo "${TEMP_LONG_CONTIG_ARRAY[*]}" )
    echo -e "\nsbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_blast_contigs.sh \
        $SCRATCH_DIR $TOOLS_DIR $BLAST_DB_TYPES $NCBI_DB_DIR_PREFIX $IDENTIFY $NUM_ALIGNMENTS $TEMP_LONG_CONTIGS_STRING\n" >> $PIPELINE_STATUS
    DEPENDENCY=$(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_blast_contigs.sh \
        $SCRATCH_DIR $TOOLS_DIR $BLAST_DB_TYPES $NCBI_DB_DIR_PREFIX $IDENTIFY $NUM_ALIGNMENTS $TEMP_LONG_CONTIGS_STRING)
    TEMP_ARRAY_START=$(($TEMP_ARRAY_START + $TEMP_ARRAY_INCREMENT))
    echo -e "$(date)\nNew start: $TEMP_ARRAY_START\nIncrement: $TEMP_ARRAY_INCREMENT" >> $PIPELINE_STATUS

    if [ $TEMP_ARRAY_START -le ${#FASTQ_ARRAY[@]} ]; then
        echo -e "\nsbatch --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        sbatch --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCY[*]}" ) -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}
        echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
    else
        echo -e "\nsbatch --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step3 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        sbatch --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step3 ${OPTIONS[@]}
        echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
    fi
elif [ $STEP -eq 3 ]; then
    BLAST_DB_TYPES_ARRAY=( $(echo $BLAST_DB_TYPES | sed 's/-/ /g') )
    BLAST_RESULTS_COUNT=$(ls ${IDENTIFY}_blast_results_* | wc -l)
    CONTIG_NUM_ARRAY=( $(ls ${IDENTIFY}_long_contigs_* | sed "s/${IDENTIFY}_long_contigs_//" | sed "s/.fasta//") )
    MAX_RESULTS=$(echo ${#CONTIG_NUM_ARRAY[@]} ${#BLAST_DB_TYPES_ARRAY[@]} | awk '{ print $1 * $2 }')
    CONTIG_COUNT=1
    for CONTIG_NUM in ${CONTIG_NUM_ARRAY[@]}; do
        for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
            if [ ! -f ${IDENTIFY}_blast_results_${DB_TYPE}_${CONTIG_NUM}.json ]; then
                echo -e "\tLong contig file number $CONTIG_COUNT - ${IDENTIFY}_blast_results_${DB_TYPE}_${CONTIG_NUM}.json not found" >> $PIPELINE_STATUS
            else
                rm ${STD_ERR_OUT_DIR}/*_${CONTIG_COUNT}_2_blast_contigs.out ${STD_ERR_OUT_DIR}/*_${CONTIG_COUNT}_2_blast_contigs.err
            fi
        done
        CONTIG_COUNT=$((CONTIG_COUNT+1))
    done
    if [ $BLAST_RESULTS_COUNT -eq 0 ]; then
        echo "No BLAST results found. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    else
        echo "$BLAST_RESULTS_COUNT BLAST results out of a possible $MAX_RESULTS maximum" >> $PIPELINE_STATUS
    fi
    echo "### BLAST aligning contigs ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    ml python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36
    echo "### Parsing BLAST results ### - START: $(date)" >> $PIPELINE_STATUS
    BLAST_DB_TYPES_ARRAY=( $(echo $BLAST_DB_TYPES | sed 's/-/ /g') )
    for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
        echo -e "\npython3 ${SCRIPT_DIR}/parse_blast_results.py $DB_TYPE \
            ${IDENTIFY}_blast_results_${DB_TYPE}_ ${IDENTIFY}.${DB_TYPE}.blast_results.tsv $PIPELINE_STATUS\n" >> $PIPELINE_STATUS
        python3 ${SCRIPT_DIR}/parse_blast_results.py $DB_TYPE \
            ${IDENTIFY}_blast_results_${DB_TYPE}_ ${IDENTIFY}.${DB_TYPE}.blast_results.tsv $PIPELINE_STATUS
    done
    echo "### Parsing BLAST results ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    echo "### Converting Kraken reports to TSV ### - START: $(date)" >> $PIPELINE_STATUS
    KRAKEN_DB_TYPE_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
    echo "Samples string: $SAMPLES_STRING"
    for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
        echo -e "\npython3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k $DB_TYPE\n" >> $PIPELINE_STATUS
        python3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k $DB_TYPE
    done
    if [ $ADD_GENUS -eq 1 ]; then
        echo -e "\npython3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k microbial\n" >> $PIPELINE_STATUS
        python3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k "microbial"
    fi
    echo "### Converting Kraken reports to TSV ### - END: $(date)" >> $PIPELINE_STATUS
    

    ml R/4.2.0
    export R_LIBS="/home/groups/cgawad/R_LIBS"
    echo "### Processing contamination, Kraken results, BLAST results, and making final figures ### - START: $(date)" >> $PIPELINE_STATUS
    echo -e "\nRscript ${SCRIPT_DIR}/analyze_and_plot_results.R \
        --project $PROJECT --identify $IDENTIFY \
        --sample_read_count_filename ${PROJECT}.sample_read_counts.tsv \
        --summed_read_targets_filename ${PROJECT}.summed_read_targets.tsv \
        --contig_read_targets_filename ${PROJECT}.contig_read_targets.tsv \
        --contig_data_filename ${PROJECT}.contig_data.tsv \
        --kraken_db_types $KRAKEN_DB_TYPES --kraken_jtree_suffix .kraken_jtree.json \
        --blast_db_types $BLAST_DB_TYPES --blast_results_suffix .blast_results.tsv \
        --contig_alignment_fraction_min $CONTIG_ALIGN_MINIMUM \
        --ncbi_annotations_dir $NCBI_ANNOTATIONS_DIR ${FIGURE_OPTIONS[@]}\n" >> $PIPELINE_STATUS
    Rscript ${SCRIPT_DIR}/analyze_and_plot_results.R \
        --project $PROJECT --identify $IDENTIFY \
        --sample_read_count_filename ${PROJECT}.sample_read_counts.tsv \
        --summed_read_targets_filename ${PROJECT}.summed_read_targets.tsv \
        --contig_read_targets_filename ${PROJECT}.contig_read_targets.tsv \
        --contig_data_filename ${PROJECT}.contig_data.tsv \
        --kraken_db_types $KRAKEN_DB_TYPES --kraken_jtree_suffix ".kraken_jtree.json" \
        --blast_db_types $BLAST_DB_TYPES --blast_results_suffix ".blast_results.tsv" \
        --contig_alignment_fraction_min $CONTIG_ALIGN_MINIMUM \
        --ncbi_annotations_dir $NCBI_ANNOTATIONS_DIR ${FIGURE_OPTIONS[@]}
    echo "### Processing contamination, Kraken results, BLAST results, and making final figures ### - END: $(date)" >> $PIPELINE_STATUS
    
    echo "### Removing intermediate files ### - START: $(date)" >> $PIPELINE_STATUS
    rm ${IDENTIFY}_long_contigs_*
    rm ${IDENTIFY}.*.kraken_jtree.json 
    rm ${IDENTIFY}_blast_results_*.json
    echo "### Removing intermediate files ### - END: $(date)" >> $PIPELINE_STATUS
    
    if [ "$SCRATCH_DIR" != "$RESULTS_DIR" ]; then
        echo "### Moving results from scratch dir to results dir ### - START: $(date)"
        rsync -ar $SCRATCH_DIR/ $RESULTS_DIR/
        echo "### Moving results from scratch dir to results dir ### - END: $(date)"
    fi
    echo "END: $(date)" >> $PIPELINE_STATUS
fi
