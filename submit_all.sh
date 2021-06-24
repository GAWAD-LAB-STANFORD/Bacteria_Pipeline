#!/bin/bash
#
#SBATCH --job-name=submit_all
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=5:00:00
#SBATCH --partition=cgawad

PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
HELP="\
Purpose: \n\t\
    This pipeline is built to identify bacterial species from pair-end fastq.gz files and remove human contamination \n\n\
Required arguments: -p/--project <arg> and either -f/--fastq_dir <arg> or -r/--results_dir <arg> \n\
Optional arguments: -b/--run_dir <arg>, --sample_sheet <arg>, --R1_suffix <arg>, --R2_suffix <arg>, --err_out_dir <arg>, \n\t\
    --skip_trimming, --rna, --skip_identify, --only_identify, --contig_len_min <arg>, --contig_align_min <arg>, --slurm <arg> \n\
Defaults: \n\t\
    If no fastq_dir specified, uses results_dir \n\t\
    If no results_dir specified, makes new directory in fastq_dir \n\t\
    sample_sheet: SampleSheet.csv \n\t\
    R1_suffix: _L001_R1_001.fastq.gz or _R1_001.fastq.gz \n\t\
    R2_suffix: _L001_R2_001.fastq.gz or _R1_001.fastq.gz \n\t\
    contig_len_min: 5000 \n\t\
    contig_align_min: 0.9 \n\n\
Run after demultiplexing: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --fastq_dir /oak/stanford/groups/cgawad/MRD_project/ --project MRD_project \n\n\
Run with demultiplexing: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/191126_MN01236_0003_A000H2WWHT --fastq_dir /oak/stanford/groups/cgawad/MRD_project/ --project MRD_project \n\n\
Run with demultiplexing, wait 12 hours before starting, and email notification when analysis begins and ends: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/191126_MN01236_0003_A000H2WWHT --fastq_dir /oak/stanford/groups/cgawad/MRD_project/ --project MRD_project --slurm --begin=now+12hours --mail-type=ALL \n\n\
For more information, read the README.md"

# Reads in command line option arguments and assigns them to variables
SKIP_TRIMMOMATIC=0
RNA_INPUT=0
SKIP_IDENTIFY=0
ONLY_IDENTIFY=0
CONTIG_LENGTH_MINIMUM=5000
CONTIG_ALIGN_MINIMUM=0.9
STEP=0
DEPENDENCIES=()
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
        --skip_trimming )       SKIP_TRIMMOMATIC=1
                                ;;
        --rna )                 RNA_INPUT=1
                                ;;
        --skip_identify )       SKIP_IDENTIFY=1
                                ;;
        --only_identify )       ONLY_IDENTIFY=1
                                ;;
        --identify )            shift
                                IDENTIFY=$1
                                ;;
        --contig_len_min )      shift
                                CONTIG_LENGTH_MINIMUM=$1
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
REFERENCE_DIR="/oak/stanford/groups/cgawad/Reference_Files/GATK_Resource_Bundle_hg38"
TOOLS_DIR="/oak/stanford/groups/cgawad/Sequencing_Analysis_Tools"
REF_FASTA="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
SCRIPT_DIR="${PIPELINE_DIR}/scripts"
KRAKEN_DB_TYPES="microbial-plasmid-viral"
KRAKEN_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/Kraken2_Fatfree_Databases/kraken2-fatfree-"
BLAST_DB_TYPES="nt-plasmid-viral"
NCBI_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/NCBI_RefSeq_Databases/ncbi_database_"
NCBI_ANNOTATIONS_DIR="/oak/stanford/groups/cgawad/Reference_Files/NCBI_Annotations"
CONTIGS_PER_BLAST_JOB=320

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
OPTIONS=( "--err_out_dir $STD_ERR_OUT_DIR -f $FASTQ_DIR -r $RESULTS_DIR -d $PIPELINE_DIR -p $PROJECT" )
if [ ! -z $RUN_DIR ] && [ $ONLY_IDENTIFY -eq 1 ]; then
    echo "Variables not supplied correctly. Cannot perform demultiplexing while only identifying data. Exiting with code 1"
    exit 1
fi
if [ ! -z $RUN_DIR ]; then
    SAMPLE_SHEET="${RUN_DIR}/SampleSheet.csv"
elif [ ! -z $SAMPLE_SHEET ]; then
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
if [ ! -z $R1_SUFFIX ]; then
    OPTIONS+=( "--R1_suffix $R1_SUFFIX" )
fi
if [ ! -z $R2_SUFFIX ]; then
    OPTIONS+=( "--R2_suffix $R2_SUFFIX" )
fi
if [ $SKIP_TRIMMOMATIC -eq 1 ]; then
    OPTIONS+=( "--skip_trimming" )
fi
if [ $RNA_INPUT -eq 1 ]; then
    OPTIONS+=( "--rna" )
fi
if [ $SKIP_IDENTIFY -eq 1 ] && [ $ONLY_IDENTIFY -eq 1 ]; then
    echo "Variables not supplied correctly. Please specify either --skip_identify or --only_identify, not both. Exiting with code 1"
    exit 1
elif [ $SKIP_IDENTIFY -eq 1 ]; then
    OPTIONS+=( "--skip_identify" )
elif [ $ONLY_IDENTIFY -eq 1 ]; then
    OPTIONS+=( "--only_identify" )
    if [ $STEP -eq 0 ]; then
        STEP=2
    fi
fi
if [ $STEP -eq 0 ] && [ -z $RUN_DIR ]; then
    STEP=1
fi
if [ ! -z $IDENTIFY ]; then
    OPTIONS+=( "--identify $IDENTIFY" )
else
    IDENTIFY=$PROJECT
fi
if [ $CONTIG_LENGTH_MINIMUM -ne 5000 ]; then
    OPTIONS+=( "--contig_len_min $CONTIG_LENGTH_MINIMUM" )
fi
if [ "$CONTIG_ALIGN_MINIMUM" = "0.9" ]; then
    OPTIONS+=( "--contig_align_min $CONTIG_ALIGN_MINIMUM" )
fi


TEMP_PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
if [ $ONLY_IDENTIFY -eq 1 ]; then
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${IDENTIFY}_pipeline_status.txt
else
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${PROJECT}_pipeline_status.txt
fi
cd $RESULTS_DIR
if [ "$TEMP_PIPELINE_DIR" = "$PIPELINE_DIR" ]; then
    echo -e "\nSTART: $(date)\nBacteria Pipeline\nErr out dir: $STD_ERR_OUT_DIR\nResults dir: $RESULTS_DIR\nProject: $PROJECT" >> $PIPELINE_STATUS
    # Optional variable definitions
    if [ $SKIP_TRIMMOMATIC -eq 1 ]; then
        echo "Option: Skip trimming - will not run trimmomatic" >> $PIPELINE_STATUS
    fi
    if [ $RNA_INPUT -eq 1 ]; then
        echo "Option: Expecting RNA input and will align using STAR instead of BWA" >> $PIPELINE_STATUS
    fi
    if [ $CONTIG_LENGTH_MINIMUM -eq 5000 ]; then
        echo "Default: Contig length minimum: 5000" >> $PIPELINE_STATUS
    else
        echo "Option: Contig length minimum: $CONTIG_LENGTH_MINIMUM" >> $PIPELINE_STATUS
    fi
    if [ "$CONTIG_ALIGN_MINIMUM" = "0.9" ]; then
        echo "Default: Contig align minimum: 0.9" >> $PIPELINE_STATUS
    else
        echo "Option: Contig align minimum: $CONTIG_ALIGN_MINIMUM" >> $PIPELINE_STATUS
    fi
    if [ $SKIP_IDENTIFY -eq 1 ]; then
        echo "Option: Skip identification of data - will only process the fastqs, build the contigs, and run Kraken2" >> $PIPELINE_STATUS
    fi
    if [ $ONLY_IDENTIFY -eq 1 ]; then
        echo "Option: Only identification of data - will only BLAST and filter from already built contigs" >> $PIPELINE_STATUS
    fi
    echo " " >> $PIPELINE_STATUS
fi


if [ ! -z $SLURM_OPTIONS ] || ([ $ONLY_IDENTIFY -eq 1 ] && [ "$TEMP_PIPELINE_DIR" = "$PIPELINE_DIR" ]); then
    echo "Option: Slurm - entire pipeline run will be queued with user parameters" >> $PIPELINE_STATUS
    sbatch -J $PROJECT ${SLURM_OPTIONS[@]} \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh ${OPTIONS[@]}
    exit 0
fi


if [ $STEP -ne 0 ] && [ $ONLY_IDENTIFY -eq 0 ]; then
    if [ -z $R1_SUFFIX ] || [ -z $R2_SUFFIX ]; then
        R1_SUFFIX="_L001_R1_001.fastq.gz"
        R2_SUFFIX="_L001_R2_001.fastq.gz"
        if [ $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" | wc -l) -eq 0 ]; then
            R1_SUFFIX="_R1_001.fastq.gz"
            R2_SUFFIX="_R2_001.fastq.gz"
        fi
    fi
    SAMPLE_ARRAY=( $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" -exec basename {} \; | \
        grep -v "Undetermined" | sed "s/${R1_SUFFIX}//") )
    if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
        echo "No fastq.gz files found in the fastq directory. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    fi
    if [ $STEP -eq 1 ]; then
        echo -e "Number of samples: ${#SAMPLE_ARRAY[@]}\nSamples: ${SAMPLE_ARRAY[@]}" >> $PIPELINE_STATUS
    fi
fi


if [ $STEP -eq 0 ]; then
    echo "### Demultiplexing ### - START: $(date)" >> $PIPELINE_STATUS
    echo -e "Run dir: $RUN_DIR\nSample sheet: $SAMPLE_SHEET" >> $PIPELINE_STATUS
    DEPENDENCIES+=( $(sbatch --parsable -e ${STD_ERR_OUT_DIR}/%A_%x.err -o ${STD_ERR_OUT_DIR}/%A_%x.out \
        ${SCRIPT_DIR}/0_demultiplexer.sh --run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET --fastq_dir $FASTQ_DIR \
        --pipeline_status $PIPELINE_STATUS) )
    sbatch --dependency=afterok:${DEPENDENCIES[0]} -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step1 ${OPTIONS[@]}
elif [ $STEP -eq 1 ]; then
    echo "### De novo assembling contigs and detecting contamination ### - START: $(date)" >> $PIPELINE_STATUS
    JOB_COUNT=${#SAMPLE_ARRAY[@]}
    echo "Process sample jobs to run: $JOB_COUNT" >> $PIPELINE_STATUS
    TEMP_ARRAY_INCREMENT=1000
    TEMP_ARRAY_START=1
    while [ $TEMP_ARRAY_START -le ${#SAMPLE_ARRAY[@]} ]; do
        TEMP_SAMPLE_ARRAY=( ${SAMPLE_ARRAY[@]:$(($TEMP_ARRAY_START - 1)):$TEMP_ARRAY_INCREMENT} )
        TEMP_JOB_COUNT=${#TEMP_SAMPLE_ARRAY[@]}
        echo "Submitting $TEMP_JOB_COUNT jobs for samples $TEMP_ARRAY_START to $(($TEMP_ARRAY_START + ${#TEMP_SAMPLE_ARRAY[@]} - 1))" >> $PIPELINE_STATUS
        TEMP_SAMPLES_STRING=$( IFS=$':'; echo "${TEMP_SAMPLE_ARRAY[*]}" )
        DEPENDENCIES+=( $(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
            --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/1_process_sample.sh \
            $FASTQ_DIR $RESULTS_DIR $R1_SUFFIX $R2_SUFFIX $SKIP_TRIMMOMATIC $REF_FASTA \
            $KRAKEN_DB_TYPES $KRAKEN_DB_DIR_PREFIX $TOOLS_DIR $RNA_INPUT $TEMP_SAMPLES_STRING) )
        TEMP_ARRAY_START=$(($TEMP_ARRAY_START + $TEMP_ARRAY_INCREMENT))
    done
    sbatch --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCIES[*]}" ) -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step2 ${OPTIONS[@]}
elif [ $STEP -eq 2 ] && [ $ONLY_IDENTIFY -eq 0 ]; then
    SAMPLE_COUNT=1
    for SAMPLE in ${SAMPLE_ARRAY[@]}; do
        if [ ! -f ${SAMPLE}_contigs.fasta ]; then
            echo -e "\tSample number $SAMPLE_COUNT - ${SAMPLE}_contigs.fasta file not found" >> $PIPELINE_STATUS
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
        rm ${STD_ERR_OUT_DIR}/*1_process_sample.out ${STD_ERR_OUT_DIR}/*1_process_sample.err
    fi
    echo "### De novo assembling contigs and detecting contamination ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    ml R/4.0.2
    export R_LIBS="/home/groups/cgawad/R_libs"
    bash ${SCRIPT_DIR}/summarize_metrics.sh $PIPELINE_STATUS $PROJECT $KRAKEN_DB_TYPES $RUN_DIR $SAMPLE_SHEET
    
    
    if [ $SKIP_IDENTIFY -eq 1 ]; then
        echo "Ending without identfication of data" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 0
    fi
fi


if [ $STEP -eq 2 ] || [ $STEP -eq 3 ]; then
    SAMPLE_ARRAY=( $(ls *_contigs.fasta | sed "s/_contigs.fasta//") )
    if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
        echo "No contig fasta files found in the results directory. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    fi
fi


if [ $STEP -eq 2 ]; then
    ml python/3.6.1 biology py-biopython/1.70_py27
    ml python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36


    echo "### Organzing contigs ### - START: $(date)" >> $PIPELINE_STATUS
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
    echo "BLAST jobs to run: $JOB_COUNT" >> $PIPELINE_STATUS
    TEMP_ARRAY_INCREMENT=1000
    TEMP_ARRAY_START=1
    while [ $TEMP_ARRAY_START -le ${#LONG_CONTIG_ARRAY[@]} ]; do
        TEMP_LONG_CONTIG_ARRAY=( ${LONG_CONTIG_ARRAY[@]:$(($TEMP_ARRAY_START - 1)):$TEMP_ARRAY_INCREMENT} )
        TEMP_JOB_COUNT=${#TEMP_LONG_CONTIG_ARRAY[@]}
        echo "Submitting $TEMP_JOB_COUNT jobs for samples $TEMP_ARRAY_START to $(($TEMP_ARRAY_START + ${#TEMP_LONG_CONTIG_ARRAY[@]} - 1))" >> $PIPELINE_STATUS
        TEMP_LONG_CONTIGS_STRING=$( IFS=$':'; echo "${TEMP_LONG_CONTIG_ARRAY[*]}" )
        echo "sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
            --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_blast_contigs.sh \
            $RESULTS_DIR $TOOLS_DIR $BLAST_DB_TYPES $NCBI_DB_DIR_PREFIX $IDENTIFY $TEMP_LONG_CONTIGS_STRING"
        DEPENDENCIES+=( $(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
            --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_blast_contigs.sh \
            $RESULTS_DIR $TOOLS_DIR $BLAST_DB_TYPES $NCBI_DB_DIR_PREFIX $IDENTIFY $TEMP_LONG_CONTIGS_STRING) )
        TEMP_ARRAY_START=$(($TEMP_ARRAY_START + $TEMP_ARRAY_INCREMENT))
    done
    sbatch --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCIES[*]}" ) -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step3 ${OPTIONS[@]}
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
        rm ${STD_ERR_OUT_DIR}/*2_blast_contigs.out ${STD_ERR_OUT_DIR}/*2_blast_contigs.err
    fi
    echo "### BLAST aligning contigs ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    ml python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36
    echo "### Parsing BLAST results ### - START: $(date)" >> $PIPELINE_STATUS
    BLAST_DB_TYPES_ARRAY=( $(echo $BLAST_DB_TYPES | sed 's/-/ /g') )
    for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
        python3 ${SCRIPT_DIR}/parse_blast_results.py $DB_TYPE \
            ${IDENTIFY}_blast_results_${DB_TYPE}_ ${IDENTIFY}.${DB_TYPE}.blast_results.tsv $PIPELINE_STATUS
    done
    echo "### Parsing BLAST results ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    echo "### Converting Kraken reports to TSV ### - START: $(date)" >> $PIPELINE_STATUS
    KRAKEN_DB_TYPE_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
    echo "Samples string: $SAMPLES_STRING"
    for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
        python3 ${SCRIPT_DIR}/kraken_report_to_jtree.py -p $PROJECT -k $DB_TYPE
    done
    echo "### Converting Kraken reports to TSV ### - END: $(date)" >> $PIPELINE_STATUS
    

    ml R/4.0.2
    export R_LIBS="/home/groups/cgawad/R_libs"
    echo "### Processing contamination, Kraken results, BLAST results, and making final figures ### - START: $(date)" >> $PIPELINE_STATUS
    Rscript ${SCRIPT_DIR}/analyze_and_plot_results.R \
        --project $PROJECT --identify $IDENTIFY \
        --sample_read_count_filename ${PROJECT}.sample_read_counts.tsv \
        --summed_read_targets_filename ${PROJECT}.summed_read_targets.tsv \
        --contig_read_targets_filename ${PROJECT}.contig_read_targets.tsv \
        --contig_data_filename ${PROJECT}.contig_data.tsv \
        --kraken_db_types $KRAKEN_DB_TYPES --kraken_jtree_suffix ".kraken_jtree.json" \
        --blast_db_types $BLAST_DB_TYPES --blast_results_suffix ".blast_results.tsv" \
        --contig_alignment_fraction_min $CONTIG_ALIGN_MINIMUM \
        --ncbi_annotations_dir $NCBI_ANNOTATIONS_DIR
    echo "### Processing contamination, Kraken results, BLAST results, and making final figures ### - END: $(date)" >> $PIPELINE_STATUS
    
    rm ${IDENTIFY}_long_contigs_
    rm ${IDENTIFY}.*.kraken_jtree.json 
    rm ${IDENTIFY}_blast_results_*.json
    echo "END: $(date)" >> $PIPELINE_STATUS
fi
