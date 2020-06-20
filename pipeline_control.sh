#!/bin/bash
#
#SBATCH --job-name=pipeline_control
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=2
#SBATCH --time=6-00:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
HELP="\
Purpose: \n\t\
    This pipeline is built to identify bacterial species from pair-end fastq.gz files and remove human contamination \n\n\
Get detailed help from the submit_all.sh script: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --help \n\n\
For more information, read the README.md"
PIPELINE_DIR=


# Reads in command line option arguments and assigns them to variables
CONTIG_LENGTH_MINIMUM=5000
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
        -d | --pipeline_dir )   shift
                                PIPELINE_DIR=$1
                                ;;
    esac
    shift
done
REFERENCE_DIR="/oak/stanford/groups/cgawad/Reference_Files/GATK_Resource_Bundle_hg38"
TOOLS_DIR="/oak/stanford/groups/cgawad/Sequencing_Analysis_Tools"
REF_FASTA="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
SCRIPT_DIR="${PIPELINE_DIR}/scripts"
# Tools and databases
KRAKEN_DB_TYPES="microbial-plasmid-viral"
KRAKEN_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/Kraken2_Fatfree_Databases/kraken2-fatfree-"
BLAST_DB_TYPES="nt-plasmid-viral"
NCBI_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/NCBI_RefSeq_Databases/ncbi_database_"
NCBI_ANNOTATIONS_DIR="/oak/stanford/groups/cgawad/Reference_Files/NCBI_Annotations"
# Variables
CONTIGS_PER_BLAST_JOB=320


# Ensure variables are set, set undefined variables, and enter results directory
if [ -z $FASTQ_DIR ] || [ -z $RESULTS_DIR ] || [ -z $PROJECT ] || [ -z $PIPELINE_DIR ]; then
    echo "Variables not supplied correctly. Use -h/--help options for assistance. Exiting with code 1"
    exit 1
fi
echo -e "START: $(date)\nBacteria Pipeline\nErr out dir: $STD_ERR_OUT_DIR\nFastq dir: $FASTQ_DIR\nResults dir: $RESULTS_DIR\nProject: $PROJECT\nContig len min: $CONTIG_LENGTH_MINIMUM"
cd $RESULTS_DIR

ml perl R/3.6.1

# Demultiplexing option
if [ ! -z $RUN_DIR ]; then
    echo "### Demultiplexing ### - START: $(date)"
    if [ -z $SAMPLE_SHEET ]; then
        SAMPLE_SHEET="${RUN_DIR}/SampleSheet.csv"
    fi
    if [ ! -f $SAMPLE_SHEET ]; then
        echo "Sample sheet $SAMPLE_SHEET not found. Exiting with code 1"
        exit 1
    fi
    echo -e "Run dir: $RUN_DIR\nSample sheet: $SAMPLE_SHEET"
    sbatch --wait -e $STD_ERR_OUT_DIR/%A_%x.err -o $STD_ERR_OUT_DIR/%A_%x.out \
        ${SCRIPT_DIR}/0_demultiplexer.sh --run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET --fastq_dir $FASTQ_DIR
    if [ $(find $FASTQ_DIR -name "*.fastq.gz" | wc -l) -eq 0 ]; then
        echo "No fastq.gz files found in fastq dir. Exiting with code 1"
        exit 1
    fi
    find $FASTQ_DIR -name "*.fastq.gz" -exec mv {} ${FASTQ_DIR}/ \;
    BCL_SIZE=$(du -sh $RUN_DIR | cut -f 1)
    UND_SIZE=$(du -shc ${FASTQ_DIR}/Undetermined*.fastq.gz | tail -n 1 | cut -f 1)
    FASTQ_SIZE=$(ls ${FASTQ_DIR}/*.fastq.gz | grep -v "Undetermined\|extracted" | xargs du -shc | tail -n 1 | cut -f 1)
    echo -e "$BCL_SIZE run dir produced $UND_SIZE of undetermined fastq.gz and $FASTQ_SIZE of determined fastq.gz"
    echo "### Demultiplexing ### - END: $(date)"
fi


if [ -z $R1_SUFFIX ] || [ -z $R2_SUFFIX ]; then
    R1_SUFFIX="_L001_R1_001.fastq.gz"
    R2_SUFFIX="_L001_R2_001.fastq.gz"
    if [ $(find ${FASTQ_DIR}/ -maxdepth 1 -name "*${R1_SUFFIX}" | wc -l) -eq 0 ]; then
        R1_SUFFIX="_R1_001.fastq.gz"
        R2_SUFFIX="_R2_001.fastq.gz"
    fi
fi
SAMPLE_ARRAY=( $(find ${FASTQ_DIR}/ -maxdepth 1 -name "*${R1_SUFFIX}" -exec basename {} \; | \
    grep -v "Undetermined" | sed "s/${R1_SUFFIX}//") )
if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
    echo "No fastq.gz files found in the sample directory. Exiting with code 1"
    exit 1
else
    echo -e "Number of samples: ${#SAMPLE_ARRAY[@]}\nSamples: ${SAMPLE_ARRAY[@]}"
fi


echo "### De novo assembling contigs and detecting contamination ### - START: $(date)"
JOB_COUNT=${#SAMPLE_ARRAY[@]}
SAMPLES_STRING=$( IFS=$':'; echo "${SAMPLE_ARRAY[*]}" )
echo -e "Assemble jobs to run: $JOB_COUNT"
sbatch --wait -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
    --array=1-${JOB_COUNT} ${SCRIPT_DIR}/1_process_sample.sh \
    $FASTQ_DIR $RESULTS_DIR $R1_SUFFIX $R2_SUFFIX $REF_FASTA \
    $KRAKEN_DB_TYPES $KRAKEN_DB_DIR_PREFIX $TOOLS_DIR $SAMPLES_STRING
if [ $(ls *_contigs.fasta | wc -l) -eq 0 ]; then
    echo "No contig files found. Exiting with code 1"
    exit 1
fi
echo "### De novo assembling contigs and detecting contamination ### - END: $(date)"


ml python/3.6.1 biology py-biopython/1.70_py27
ml python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36


echo "### Summarizing metrics ### - START: $(date)"
SAMPLE_READ_COUNTS="${PROJECT}.sample_read_counts.tsv"
READ_COUNT_FILENAMES=( $(ls *.read_counts.tsv) )
head -n 1 ${READ_COUNT_FILENAMES[0]} > $SAMPLE_READ_COUNTS
for i in ${READ_COUNT_FILENAMES[@]}; do tail -n +2 $i; done >> $SAMPLE_READ_COUNTS
if [ ! -z $RUN_DIR ]; then
    DESIRED_CLUSTER_DENSITY=$(sed -n $(grep -n "Lane" $SAMPLE_SHEET | sed "s/:Lane.*//" | \
        awk '{print $1 + 1}')p $SAMPLE_SHEET | cut -d , -f $(grep "Lane" $SAMPLE_SHEET | \
        sed "s/,/\n/g" | nl | grep "Desired_Cluster_Density" | cut -f 1))
    INITAL_CONCENTRATION_COLUMN=$(grep "Initial_Concentration" $SAMPLE_SHEET | tr , "\n" | nl | grep "Initial_Concentration" | cut -f 1)
    LIBRARY_GROUP_COLUMN=$(grep "Library_Group" $SAMPLE_SHEET | tr , "\n" | nl | grep "Library_Group" | cut -f 1)
    if [ ! -z $INITAL_CONCENTRATION_COLUMN ] && [ ! -z $DESIRED_CLUSTER_DENSITY ]; then
        echo "### Calculating library concentration corrections ### - START: $(date)"
        cat $SAMPLE_READ_COUNTS > ${SAMPLE_READ_COUNTS}.temp
        tail -n +$(cat -n $SAMPLE_SHEET | grep "Initial_Concentration" | cut -f 1 | tr -d '[:blank:]') $SAMPLE_SHEET | \
            cut -d , -f $INITAL_CONCENTRATION_COLUMN | awk '{print tolower($0)}' | \
            paste ${SAMPLE_READ_COUNTS}.temp - > $SAMPLE_READ_COUNTS
        if [ ! -z $LIBRARY_GROUP_COLUMN ]; then
            cat $SAMPLE_READ_COUNTS > ${SAMPLE_READ_COUNTS}.temp
            tail -n +$(cat -n $SAMPLE_SHEET | grep "Library_Group" | cut -f 1 | tr -d '[:blank:]') $SAMPLE_SHEET | \
                cut -d , -f $LIBRARY_GROUP_COLUMN | awk '{print tolower($0)}' | \
                paste ${SAMPLE_READ_COUNTS}.temp - > $SAMPLE_READ_COUNTS
        fi
        rm ${SAMPLE_READ_COUNTS}.temp
        Rscript ${SCRIPT_DIR}/correct_library_concentrations.R \
            "${RUN_DIR}/RunCompletionStatus.xml" $DESIRED_CLUSTER_DENSITY $SAMPLE_READ_COUNTS $PROJECT
        echo "### Calculating library concentration corrections ### - END: $(date)"
    fi
fi

HUMAN_ALIGNMENT_METRICS_FILENAMES=( $(ls *_human_alignment_metrics.tsv) )
echo -e sample"\t"$(head -n 7 ${HUMAN_ALIGNMENT_METRICS_FILENAMES[0]} | tail -n 1) | sed 's/ /\t/g' > ${PROJECT}.human_alignment_metrics.tsv
for i in ${HUMAN_ALIGNMENT_METRICS_FILENAMES[@]}; do 
    SAMPLE=$(echo $i | sed "s/_human_alignment_metrics.tsv//")
    R1=$(head -n 8 $i | tail -n 1)
    R2=$(head -n 9 $i | tail -n 1)
    PAIR=$(head -n 10 $i | tail -n 1)
    echo -e "$SAMPLE\t$R1\n$SAMPLE\t$R2\n$SAMPLE\t$PAIR"
done | sed 's/ /\t/g' >> ${PROJECT}.human_alignment_metrics.tsv
echo "Merged human alignment metrics"

CONTIG_ALIGNMENT_METRICS_FILENAMES=( $(ls *_contig_alignment_metrics.tsv) )
echo -e sample"\t"$(head -n 7 ${CONTIG_ALIGNMENT_METRICS_FILENAMES[0]} | tail -n 1) | sed 's/ /\t/g' > ${PROJECT}.contig_alignment_metrics.tsv
for i in ${CONTIG_ALIGNMENT_METRICS_FILENAMES[@]}; do 
    SAMPLE=$(echo $i | sed "s/_contig_alignment_metrics.tsv//")
    R1=$(head -n 8 $i | tail -n 1)
    R2=$(head -n 9 $i | tail -n 1)
    PAIR=$(head -n 10 $i | tail -n 1)
    echo -e "$SAMPLE\t$R1\n$SAMPLE\t$R2\n$SAMPLE\t$PAIR"
done | sed 's/ /\t/g' >> ${PROJECT}.contig_alignment_metrics.tsv
echo "Merged contig alignment metrics"

KRAKEN_DB_TYPES_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
for DB_TYPE in ${KRAKEN_DB_TYPES_ARRAY[@]}; do
    echo -e "sample\tpercent_fragments_covered\tfragments_covered\tfragments_assigned\trank_code\ttaxid\tsciname" > \
        ${PROJECT}.${DB_TYPE}.kraken_reports.tsv
    KRAKEN_REPORT_FILENAMES=( $(ls *_${DB_TYPE}_kraken_report.tsv) )
    for i in ${KRAKEN_REPORT_FILENAMES[@]}; do
        SAMPLE=$(echo $i | sed "s/_${DB_TYPE}_kraken_report.tsv//")
        cat $i | sed 's/^ \+/'${SAMPLE}'\t/' | awk '$2>=1.00' >> ${PROJECT}.${DB_TYPE}.kraken_reports.tsv
    done
    python3 ${SCRIPT_DIR}/kraken_report_to_jtree.py ${PROJECT}.${DB_TYPE}.kraken_reports.tsv \
        $PROJECT .${DB_TYPE}.kraken_jtree.json
done
KRAKEN_REPORT_FILENAMES=( $(ls *_kraken_report.tsv) )
echo "Merged kraken reports"

echo -e "sample\tfasta_header" > ${PROJECT}.contig_data.tsv
CONTIGS_FILENAMES=( $(ls *_contigs.fasta) )
for i in ${CONTIGS_FILENAMES[@]}; do 
    grep ">" $i | xargs -i echo -e $(echo $i | sed "s/_contigs.fasta//")"\t"{} >> ${PROJECT}.contig_data.tsv; 
done
echo "Merged contigs"

SUMMED_READ_TARGETS_FILENAMES=( $(ls *_summed_read_targets.tsv) )
head -n 1 ${SUMMED_READ_TARGETS_FILENAMES[0]} > ${PROJECT}.summed_read_targets.tsv
for i in ${SUMMED_READ_TARGETS_FILENAMES[@]}; do sed -n 2p $i >> ${PROJECT}.summed_read_targets.tsv; done
echo "Merged summarized read targets"

CONTIG_READ_TARGETS_FILENAMES=( $(ls *_contig_read_targets.tsv) )
head -n 1 ${CONTIG_READ_TARGETS_FILENAMES[0]} > ${PROJECT}.contig_read_targets.tsv
for i in ${CONTIG_READ_TARGETS_FILENAMES[@]}; do tail -n +2 $i >> ${PROJECT}.contig_read_targets.tsv; done
echo "Merged contig read targets"
echo "### Summarizing metrics ### - END: $(date)"


echo "### Organzing contigs ### - START: $(date)"
python3 ${SCRIPT_DIR}/organize_contigs.py \
    "_contigs.fasta" $CONTIG_LENGTH_MINIMUM $CONTIGS_PER_BLAST_JOB "long_contigs_"
if [ $(ls long_contigs_* | wc -l) -eq 0 ]; then
    echo "No contigs longer than $CONTIG_LENGTH_MINIMUM. Exiting with code 1"
    exit 1
fi
echo "### Organzing contigs ### - END: $(date)"


echo "### BLAST aligning contigs to nucleotide database ### - START: $(date)"
JOB_COUNT=$(ls long_contigs_* | wc -l)
echo -e "Blast jobs to run: $JOB_COUNT"
sbatch --wait -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
    --array=1-${JOB_COUNT} ${SCRIPT_DIR}/2_blast_contigs.sh \
    $RESULTS_DIR $TOOLS_DIR $BLAST_DB_TYPES $NCBI_DB_DIR_PREFIX
if [ $(ls blast_results_* | wc -l) -eq 0 ]; then
    echo "No BLAST results found. Exiting with code 1"
    exit 1
fi
echo "### BLAST aligning contigs to nucleotide database ### - END: $(date)"


echo "### Parsing BLAST results ### - START: $(date)"
BLAST_DB_TYPES_ARRAY=( $(echo $BLAST_DB_TYPES | sed 's/-/ /g') )
for DB_TYPE in ${BLAST_DB_TYPES_ARRAY[@]}; do
    python3 ${SCRIPT_DIR}/parse_blast_results.py $DB_TYPE \
        blast_results_${DB_TYPE}_ ${PROJECT}.${DB_TYPE}.blast_results.tsv
done
echo "### Parsing BLAST results ### - END: $(date)"


echo "### Processing contamination and BLAST results ### - START: $(date)"
Rscript ${SCRIPT_DIR}/analyze_and_plot_results.R \
    $PROJECT $SAMPLE_READ_COUNTS ${PROJECT}.summed_read_targets.tsv \
    ${PROJECT}.contig_read_targets.tsv ${PROJECT}.contig_data.tsv \
    $KRAKEN_DB_TYPES ".kraken_jtree.json" \
    $BLAST_DB_TYPES ".blast_results.tsv" $NCBI_ANNOTATIONS_DIR
echo "### Processing contamination and BLAST results ### - END: $(date)"


if [ ! -f ${PROJECT}.${BLAST_DB_TYPES_ARRAY[0]}.blast_results.tsv ]; then
    echo "Final file ${PROJECT}.${BLAST_DB_TYPES_ARRAY[0]}.blast_results.tsv not found. Exiting with code 1"
    exit 1
fi
rm ${READ_COUNT_FILENAMES[@]} ${PROJECT}.*.kraken_jtree.json
rm ${KRAKEN_REPORT_FILENAMES[@]} ${CONTIGS_FILENAMES[@]} 
rm ${SUMMED_READ_TARGETS_FILENAMES[@]} ${CONTIG_READ_TARGETS_FILENAMES[@]}
rm ${HUMAN_ALIGNMENT_METRICS_FILENAMES[@]} ${CONTIG_ALIGNMENT_METRICS_FILENAMES[@]}
rm blast_results_*.json # long_contigs_*
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"