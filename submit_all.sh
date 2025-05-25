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
    --analyze_scaffolds, --query_len_min <arg>, --kraken_db_types <arg>, --blast_db_types <arg>, \n\t\
    --blast_num_alignments <arg>, --blast_align_min <arg>, --blast_hit_rank_min <arg>, --query_align_min <arg>, --add_genus, --slurm <arg> \n\
Defaults: \n\t\
    If no fastq_dir specified, uses results_dir \n\t\
    If no results_dir specified, makes new directory in fastq_dir \n\t\
    scratch_dir: /scratch/groups/cgawad/date_project_Scratch \n\t\
    sample_sheet: SampleSheet.csv \n\t\
    R1_suffix: _L001_R1_001.fastq.gz or _R1_001.fastq.gz or _R1.fastq.gz \n\t\
    R2_suffix: _L001_R2_001.fastq.gz or _R2_001.fastq.gz or _R2.fastq.gz \n\t\
    query_len_min: 5000 \n\t\
    kraken_db_types: microbial \n\t\
    blast_db_types: nt \n\t\
    blast_num_alignments: 5 \n\t\
    blast_align_min: 0 \n\t\
    blast_hit_rank_min: 1 \n\t\
    query_align_min: 0.9 \n\n\
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
QUERY="contig"
QUERY_LENGTH_MIN=5000
KRAKEN_DB_TYPES="microbial"
BLAST_DB_TYPES="nt"
BLAST_NUM_ALIGNMENTS=5
BLAST_ALIGN_MIN=0
BLAST_HIT_RANK_MIN=1
QUERY_ALIGN_MIN=0.9
ADD_GENUS=0
STEP=0
TEMP_ARRAY_START=0
DEPENDENCY=""
DEPENDER=""
FIGURE_OPTIONS=()
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )               echo -e $HELP
                                    exit 0
                                    ;;
        -f | --fastq_dir )          shift
                                    FASTQ_DIR=$1
                                    ;;
        -r | --results_dir )        shift
                                    RESULTS_DIR=$1
                                    ;;
        -p | --project )            shift
                                    PROJECT=$1
                                    ;;
        -d | --pipeline_dir )       shift
                                    PIPELINE_DIR=$1
                                    ;;
        -s | --scratch_dir )        shift
                                    SCRATCH_DIR=$1
                                    ;;
        --err_out_dir )             shift
                                    STD_ERR_OUT_DIR=$1
                                    ;;
        --skip_scratch )            SKIP_SCRATCH=1
                                    ;;
        -b | --run_dir )            shift
                                    RUN_DIR=$1
                                    ;;
        --sample_sheet )            shift
                                    SAMPLE_SHEET=$1
                                    ;;
        --R1_suffix )               shift
                                    R1_SUFFIX=$1
                                    ;;
        --R2_suffix )               shift
                                    R2_SUFFIX=$1
                                    ;;
        --skip_identify )           SKIP_IDENTIFY=1
                                    ;;
        --only_identify )           ONLY_IDENTIFY=1
                                    ;;
        --identify )                shift
                                    IDENTIFY=$1
                                    ;;
        --skip_trimming )           SKIP_TRIMMOMATIC=1
                                    ;;
        --rna )                     RNA=1
                                    ;;
        --filter_rhesus )           FILTER_RHESUS=1
                                    ;;
        --analyze_scaffolds )       QUERY="scaffold"
                                    ;;
        --query_len_min )           shift
                                    QUERY_LENGTH_MIN=$1
                                    ;;
        --add_genus )               shift
                                    ADD_GENUS=1
                                    ;;
        --kraken_db_types )         shift
                                    KRAKEN_DB_TYPES=$1
                                    ;;
        --blast_db_types )          shift
                                    BLAST_DB_TYPES=$1
                                    ;;
        --blast_num_alignments )    shift
                                    BLAST_NUM_ALIGNMENTS=$1
                                    ;;
        --blast_align_min )         shift
                                    BLAST_ALIGN_MIN=$1
                                    ;;
        --blast_hit_rank_min )      shift
                                    BLAST_HIT_RANK_MIN=$1
                                    ;;
        --query_align_min )         shift
                                    QUERY_ALIGN_MIN=$1
                                    ;;
        --step1 )                   STEP=1
                                    ;;
        --step2 )                   STEP=2
                                    ;;
        --step3 )                   STEP=3
                                    ;;
        --temp_array_start )        shift
                                    TEMP_ARRAY_START=$1
                                    ;;
        --slurm )                   shift
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
QUERIES_PER_BLAST_JOB=320

# Ensure we have the required variables set and set other variables
if ([ -z $FASTQ_DIR ] && [ -z $RESULTS_DIR ]) || [ -z $PROJECT ] || [ -z $PIPELINE_DIR ]; then
    echo "Variables not supplied correctly. Use -h/--help options for assistance. Exiting with code 1"
    exit 1
fi
if ([ -z $FASTQ_DIR ] && [ $ONLY_IDENTIFY -eq 0 ]); then
    FASTQ_DIR="$RESULTS_DIR"
    OPTIONS+=( "-f $FASTQ_DIR" )
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
mkdir -p $FASTQ_DIR
mkdir -p $RESULTS_DIR
mkdir -p $SCRATCH_DIR
mkdir -p $STD_ERR_OUT_DIR

# Add parameters/arguments to OPTIONS variable to retain them with each subsequent pipeline resubmission
OPTIONS=( "-r $RESULTS_DIR -d $PIPELINE_DIR -p $PROJECT -s $SCRATCH_DIR --err_out_dir $STD_ERR_OUT_DIR " )
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
if [ ! -z $FASTQ_DIR ]; then
    OPTIONS+=( "-f $FASTQ_DIR" )
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
if [ "$QUERY" != "contig" ]; then
    OPTIONS+=( "--analyze_scaffolds" )
fi
if [ $QUERY_LENGTH_MIN -ne 5000 ]; then
    OPTIONS+=( "--query_len_min $QUERY_LENGTH_MIN" )
fi
if [ "$KRAKEN_DB_TYPES" != "microbial" ]; then
    OPTIONS+=( "--kraken_db_types $KRAKEN_DB_TYPES" )
fi
if [ "$BLAST_DB_TYPES" != "nt" ]; then
    OPTIONS+=( "--blast_db_types $BLAST_DB_TYPES" )
fi
if [ $BLAST_NUM_ALIGNMENTS -ne 250 ]; then
    OPTIONS+=( "--blast_num_alignments $BLAST_NUM_ALIGNMENTS" )
fi
if [ "$BLAST_HIT_RANK_MIN" != "1" ]; then
    OPTIONS+=( "--blast_hit_rank_min $BLAST_HIT_RANK_MIN" )
fi
if [ "$QUERY_ALIGN_MIN" != "0.9" ]; then
    OPTIONS+=( "--query_align_min $QUERY_ALIGN_MIN" )
fi
if [ $ADD_GENUS -eq 1 ]; then
    FIGURE_OPTIONS+=( "--add_genus" )
fi

# On first run, stdout all options from OPTIONS variable for pipeline resubmission parameters/arguments
TEMP_PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
if [ $ONLY_IDENTIFY -eq 1 ]; then
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${IDENTIFY}_pipeline_status.txt
else
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${PROJECT}_pipeline_status.txt
fi
cd $SCRATCH_DIR
if [ "$TEMP_PIPELINE_DIR" = "$PIPELINE_DIR" ]; then
    echo -e "\nSTART: $(date)\nBacteria Pipeline\n\n$PIPELINE_DIR/submit_all.sh $PIPELINE_COMMAND\n\nProject: $PROJECT\nResults dir: $RESULTS_DIR\nScratch dir: $SCRATCH_DIR\nErr out dir: $STD_ERR_OUT_DIR\n" >> $PIPELINE_STATUS
    echo -e "\nOPTIONS variable holding parameters/arguments for pipeline resubmission: ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
fi


if [ ! -z $SLURM_OPTIONS ]; then
    echo "Option: Slurm - entire pipeline run will be queued with user parameters" >> $PIPELINE_STATUS
    sbatch -J $PROJECT ${SLURM_OPTIONS[@]} \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh ${OPTIONS[@]}
    exit 0
fi


if ([ $STEP -eq 0 ] && [ -z $RUN_DIR ] && [ $ONLY_IDENTIFY -eq 0 ]) || [ $STEP -eq 1 ] || [ $STEP -eq 2 ]; then
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
    echo -e "\nsbatch --parsable --dependency=afterok:$DEPENDENCY -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step1 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
    DEPENDER=$(sbatch --parsable --dependency=afterok:$DEPENDENCY -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step1 ${OPTIONS[@]})
    echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
elif ([ $STEP -eq 0 ] && [ -z $RUN_DIR ] && [ $ONLY_IDENTIFY -eq 0 ]) || [ $STEP -eq 1 ]; then
    if [ $TEMP_ARRAY_START -eq 0 ]; then
        echo -e "Number of samples: ${#SAMPLE_ARRAY[@]}\nSamples: ${SAMPLE_ARRAY[@]}\n" >> $PIPELINE_STATUS
        echo "### De novo assembling contigs and scaffolds, and detecting contamination ### - START: $(date)" >> $PIPELINE_STATUS
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
        --fastq_dir $FASTQ_DIR --scratch_dir $SCRATCH_DIR --R1_suffix $R1_SUFFIX --R2_suffix $R2_SUFFIX \
        --skip_trimmomatic $SKIP_TRIMMOMATIC --rna $RNA --ref_fasta_string $REF_FASTA_STRING \
        --ref_name_string $REF_NAME_STRING --kraken_db_types_string $KRAKEN_DB_TYPES \
        --kraken_db_dir_prefix $KRAKEN_DB_DIR_PREFIX --tools_dir $TOOLS_DIR --sample_string $TEMP_SAMPLES_STRING\n" >> $PIPELINE_STATUS
    DEPENDENCY=$(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/1_process_sample.sh \
        --fastq_dir $FASTQ_DIR --scratch_dir $SCRATCH_DIR --R1_suffix $R1_SUFFIX --R2_suffix $R2_SUFFIX \
        --skip_trimmomatic $SKIP_TRIMMOMATIC --rna $RNA --ref_fasta_string $REF_FASTA_STRING \
        --ref_name_string $REF_NAME_STRING --kraken_db_types_string $KRAKEN_DB_TYPES \
        --kraken_db_dir_prefix $KRAKEN_DB_DIR_PREFIX --tools_dir $TOOLS_DIR --sample_string $TEMP_SAMPLES_STRING)
    TEMP_ARRAY_START=$(($TEMP_ARRAY_START + $TEMP_ARRAY_INCREMENT))
    echo -e "$(date)\nIncrement: $TEMP_ARRAY_INCREMENT\nNew start: $TEMP_ARRAY_START" >> $PIPELINE_STATUS

    if [ $TEMP_ARRAY_START -le ${#SAMPLE_ARRAY[@]} ]; then
        echo -e "\nsbatch --parsable --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step1 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        DEPENDER=$(sbatch --parsable --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step1 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]})
    else
        echo -e "\nsbatch --parsable --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        DEPENDER=$(sbatch --parsable --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 ${OPTIONS[@]})
    fi
    echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
elif [ $STEP -eq 2 ] && [ $TEMP_ARRAY_START -eq 0 ]; then
    SAMPLE_COUNT=1
    for SAMPLE in ${SAMPLE_ARRAY[@]}; do
        if [ ! -f ${SAMPLE}_contigs.fasta ]; then
            echo -e "\tSample number $SAMPLE_COUNT - ${SAMPLE}_contigs.fasta file not found" >> $PIPELINE_STATUS
        else
#            rm ${STD_ERR_OUT_DIR}/*_${SAMPLE_COUNT}_1_process_sample.out ${STD_ERR_OUT_DIR}/*_${SAMPLE_COUNT}_1_process_sample.err
		echo "next" 
       fi
        SAMPLE_COUNT=$((SAMPLE_COUNT+1))
    done
    CONTIG_FILE_COUNT=$(ls *_contigs.fasta | wc -l)
    SCAFFOLD_FILE_COUNT=$(ls *_scaffolds.fasta | wc -l)
    if [ $CONTIG_FILE_COUNT -eq 0 ]; then
        echo "No contigs found. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    else
        echo "$CONTIG_FILE_COUNT contig files out of a possible ${#SAMPLE_ARRAY[@]} maximum" >> $PIPELINE_STATUS
        echo "$SCAFFOLD_FILE_COUNT scaffold files out of a possible ${#SAMPLE_ARRAY[@]} maximum" >> $PIPELINE_STATUS
    fi
    echo "### De novo assembling contigs and scaffolds, and detecting contamination ### - END: $(date)" >> $PIPELINE_STATUS

    ml R/4.2.0
    export R_LIBS="/home/groups/cgawad/R_libs"
    bash ${SCRIPT_DIR}/merge_metrics.sh $PIPELINE_STATUS $PROJECT $KRAKEN_DB_TYPES $REF_NAME_STRING $RUN_DIR $SAMPLE_SHEET

    if [ $SKIP_IDENTIFY -eq 1 ]; then
        echo "Ending without identfication of data" >> $PIPELINE_STATUS
        if [ "$SCRATCH_DIR" != "$RESULTS_DIR" ]; then
            echo "### Moving results from scratch dir to results dir ### - START: $(date)"
            rsync -ar $SCRATCH_DIR/ $RESULTS_DIR/
            echo "### Moving results from scratch dir to results dir ### - END: $(date)"
        fi
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 0
    fi
fi


if ([ $STEP -eq 0 ] && [ $ONLY_IDENTIFY -eq 1 ]) || [ $STEP -eq 2 ]; then
    SAMPLE_ARRAY=( $(ls *_${QUERY}s.fasta | sed "s/_${QUERY}s.fasta//") )
    if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
        echo "No ${QUERY} fasta files found in the results directory. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    fi
fi


if ([ $STEP -eq 0 ] && [ $ONLY_IDENTIFY -eq 1 ]) || [ $STEP -eq 2 ]; then
    if [ $TEMP_ARRAY_START -eq 0 ]; then
        ml python/3.6.1 biology py-biopython/1.79_py39

        echo "### Organzing queries ### - START: $(date)" >> $PIPELINE_STATUS
        echo -e "\npython3 ${SCRIPT_DIR}/organize_queries.py \
            _${QUERY}s.fasta $QUERY_LENGTH_MIN $QUERIES_PER_BLAST_JOB ${IDENTIFY}_long_${QUERY}s_ $PIPELINE_STATUS" >> $PIPELINE_STATUS
        python3 ${SCRIPT_DIR}/organize_queries.py \
            "_${QUERY}s.fasta" $QUERY_LENGTH_MIN $QUERIES_PER_BLAST_JOB "${IDENTIFY}_long_${QUERY}s_" $PIPELINE_STATUS
        if [ $(ls ${IDENTIFY}_long_${QUERY}s_* | wc -l) -eq 0 ]; then
            echo "No ${QUERY}s longer than $QUERY_LENGTH_MIN. Exiting with code 1" >> $PIPELINE_STATUS
            echo "END: $(date)" >> $PIPELINE_STATUS
            exit 1
        fi
        echo "### Organzing queries ### - END: $(date)" >> $PIPELINE_STATUS

        echo "### BLAST aligning queries ### - START: $(date)" >> $PIPELINE_STATUS
        LONG_QUERY_ARRAY=( $(ls ${IDENTIFY}_long_${QUERY}s_*) )
        JOB_COUNT=${#LONG_QUERY_ARRAY[@]}
        echo "Jobs to run: $JOB_COUNT" >> $PIPELINE_STATUS
        TEMP_ARRAY_START=1
    fi


    TEMP_LONG_QUERY_ARRAY=( ${LONG_QUERY_ARRAY[@]:$(($TEMP_ARRAY_START - 1)):$TEMP_ARRAY_INCREMENT} )
    TEMP_JOB_COUNT=${#TEMP_LONG_QUERY_ARRAY[@]}
    echo "Submitting $TEMP_JOB_COUNT jobs for samples $TEMP_ARRAY_START to $(($TEMP_ARRAY_START + ${#TEMP_LONG_QUERY_ARRAY[@]} - 1))" >> $PIPELINE_STATUS
    TEMP_LONG_QUERIES_STRING=$( IFS=$':'; echo "${TEMP_LONG_QUERY_ARRAY[*]}" )
    echo -e "\nsbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_blast_queries.sh \
        --scratch_dir $SCRATCH_DIR --tools_dir $TOOLS_DIR --blast_db_types_string $BLAST_DB_TYPES \
        --ncbi_db_dir_prefix $NCBI_DB_DIR_PREFIX --identify $IDENTIFY \
        --blast_num_alignments $BLAST_NUM_ALIGNMENTS --blast_align_min $BLAST_ALIGN_MIN \
        --query $QUERY --long_query_string $TEMP_LONG_QUERIES_STRING\n" >> $PIPELINE_STATUS
    DEPENDENCY=$(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_blast_queries.sh \
        --scratch_dir $SCRATCH_DIR --tools_dir $TOOLS_DIR --blast_db_types_string $BLAST_DB_TYPES \
        --ncbi_db_dir_prefix $NCBI_DB_DIR_PREFIX --identify $IDENTIFY \
        --blast_num_alignments $BLAST_NUM_ALIGNMENTS --blast_align_min $BLAST_ALIGN_MIN \
        --query $QUERY --long_query_string $TEMP_LONG_QUERIES_STRING)
    echo -e "$(date)\nNew start: $TEMP_ARRAY_START\nIncrement: $TEMP_ARRAY_INCREMENT" >> $PIPELINE_STATUS

    if [ $TEMP_ARRAY_START -le ${#FASTQ_ARRAY[@]} ]; then
        echo -e "\nsbatch --parsable --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        DEPENDER=$(sbatch --parsable --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCY[*]}" ) -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]})
    else
        echo -e "\nsbatch --parsable --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step3 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        DEPENDER=$(sbatch --parsable --dependency=afterany:$DEPENDENCY -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step3 ${OPTIONS[@]})
    fi
    echo -e "Dependency job array number: $DEPENDENCY\nDepender job number: $DEPENDER" >> $PIPELINE_STATUS
elif [ $STEP -eq 3 ]; then
    BLAST_DB_TYPES_ARRAY=( $(echo $BLAST_DB_TYPES | sed 's/-/ /g') )
    QUERY_NUM_ARRAY=( $(ls ${IDENTIFY}_long_${QUERY}s_* | sed "s/${IDENTIFY}_long_${QUERY}s_//" | sed "s/.fasta//") )
    QUERY_COUNT=1
    for QUERY_NUM in ${QUERY_NUM_ARRAY[@]}; do
        for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
            if [ ! -f ${IDENTIFY}_blast_results_${DB_TYPE}_${QUERY_NUM}.json ]; then
                echo -e "\tLong $QUERY file number $QUERY_COUNT - ${IDENTIFY}_blast_results_${DB_TYPE}_${QUERY_NUM}.json not found" >> $PIPELINE_STATUS
            else
		echo "next"
#                rm ${STD_ERR_OUT_DIR}/*_${QUERY_COUNT}_2_blast_queries.out ${STD_ERR_OUT_DIR}/*_${QUERY_COUNT}_2_blast_queries.err
            fi
        done
        QUERY_COUNT=$((QUERY_COUNT+1))
    done
    BLAST_RESULTS_COUNT=$(ls ${IDENTIFY}_blast_results_* | wc -l)
    MAX_RESULTS=$(echo ${#QUERY_NUM_ARRAY[@]} ${#BLAST_DB_TYPES_ARRAY[@]} | awk '{ print $1 * $2 }')
    if [ $BLAST_RESULTS_COUNT -eq 0 ]; then
        echo "No BLAST results found. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    else
        echo "$BLAST_RESULTS_COUNT BLAST results out of a possible $MAX_RESULTS maximum" >> $PIPELINE_STATUS
    fi
    echo "### BLAST aligning queries ### - END: $(date)" >> $PIPELINE_STATUS

    ml python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36
    echo "### Parsing BLAST results ### - START: $(date)" >> $PIPELINE_STATUS
    BLAST_DB_TYPES_ARRAY=( $(echo $BLAST_DB_TYPES | sed 's/-/ /g') )
    for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
        echo -e "\npython3 ${SCRIPT_DIR}/parse_blast_results.py $DB_TYPE \
            ${IDENTIFY}_blast_results_${DB_TYPE}_ ${IDENTIFY}.${DB_TYPE}_${QUERY}.blast_results.tsv $PIPELINE_STATUS\n" >> $PIPELINE_STATUS
        python3 ${SCRIPT_DIR}/parse_blast_results.py $DB_TYPE \
            ${IDENTIFY}_blast_results_${DB_TYPE}_ ${IDENTIFY}.${DB_TYPE}_${QUERY}.blast_results.tsv $PIPELINE_STATUS
    done
    echo "### Parsing BLAST results ### - END: $(date)" >> $PIPELINE_STATUS

    echo "### Converting Kraken reports to TSV ### - START: $(date)" >> $PIPELINE_STATUS
    KRAKEN_DB_TYPE_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
    echo "Samples string: $SAMPLES_STRING"
    for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
        echo -e "\npython3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k $DB_TYPE -q $QUERY\n" >> $PIPELINE_STATUS
        python3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k $DB_TYPE -q $QUERY
    done
    if [ $ADD_GENUS -eq 1 ]; then
        echo -e "\npython3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k microbial -q $QUERY\n" >> $PIPELINE_STATUS
        python3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k "microbial" -q $QUERY
    fi
    echo "### Converting Kraken reports to TSV ### - END: $(date)" >> $PIPELINE_STATUS

    ml R/4.2.0
    export R_LIBS="/home/groups/cgawad/R_LIBS"
    echo "### Processing contamination, Kraken results, BLAST results, and making final figures ### - START: $(date)" >> $PIPELINE_STATUS
    echo -e "\nRscript ${SCRIPT_DIR}/analyze_and_plot_results.R \
        --project $PROJECT --identify $IDENTIFY \
        --sample_read_count_filename ${PROJECT}.sample_read_counts.tsv \
        --summed_read_targets_filename ${PROJECT}.summed_read_targets.tsv \
        --query $QUERY \
        --query_read_targets_filename ${PROJECT}.${QUERY}_read_targets.tsv \
        --query_data_filename ${PROJECT}.${QUERY}_data.tsv \
        --kraken_db_types $KRAKEN_DB_TYPES --kraken_jtree_suffix _${QUERY}.kraken_jtree.json \
        --blast_db_types $BLAST_DB_TYPES --blast_results_suffix _${QUERY}.blast_results.tsv \
        --blast_hit_rank_min $BLAST_HIT_RANK_MIN \
        --query_align_min $QUERY_ALIGN_MIN \
        --ncbi_annotations_dir $NCBI_ANNOTATIONS_DIR ${FIGURE_OPTIONS[@]}\n" >> $PIPELINE_STATUS
    Rscript ${SCRIPT_DIR}/analyze_and_plot_results.R \
        --project $PROJECT --identify $IDENTIFY \
        --sample_read_count_filename ${PROJECT}.sample_read_counts.tsv \
        --summed_read_targets_filename ${PROJECT}.summed_read_targets.tsv \
        --query $QUERY \
        --query_read_targets_filename ${PROJECT}.${QUERY}_read_targets.tsv \
        --query_data_filename ${PROJECT}.${QUERY}_data.tsv \
        --kraken_db_types $KRAKEN_DB_TYPES --kraken_jtree_suffix "_${QUERY}.kraken_jtree.json" \
        --blast_db_types $BLAST_DB_TYPES --blast_results_suffix "_${QUERY}.blast_results.tsv" \
        --blast_hit_rank_min $BLAST_HIT_RANK_MIN \
        --query_align_min $QUERY_ALIGN_MIN \
        --ncbi_annotations_dir $NCBI_ANNOTATIONS_DIR ${FIGURE_OPTIONS[@]}
    echo "### Processing contamination, Kraken results, BLAST results, and making final figures ### - END: $(date)" >> $PIPELINE_STATUS

    echo "### Removing intermediate files ### - START: $(date)" >> $PIPELINE_STATUS
#    rm ${IDENTIFY}_long_${QUERY}s_*
#    rm ${IDENTIFY}.*.kraken_jtree.json
#    rm ${IDENTIFY}_blast_results_*.json
    echo "### Removing intermediate files ### - END: $(date)" >> $PIPELINE_STATUS

    if [ "$SCRATCH_DIR" != "$RESULTS_DIR" ]; then
        echo "### Moving results from scratch dir to results dir ### - START: $(date)"
        rsync -ar $SCRATCH_DIR/ $RESULTS_DIR/
        echo "### Moving results from scratch dir to results dir ### - END: $(date)"
    fi
    echo "END: $(date)" >> $PIPELINE_STATUS
fi
