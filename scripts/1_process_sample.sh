#!/bin/bash
#
#SBATCH --job-name=1_process_sample
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
FASTQ_DIR=$1
RESULTS_DIR=$2
R1_SUFFIX=$3
R2_SUFFIX=$4
REF_FASTA=$5
CONTIGS_SUFFIX=$6
SUMMED_READ_TARGETS_SUFFIX=$7
CONTIGS_READ_TARGETS_SUFFIX=$8
KRAKEN_TOOL_DIR=$9
KRAKEN_DB_TYPES_ARRAY=( $(echo ${10} | sed 's/-/ /g') )
KRAKEN_DB_DIR_PREFIX=${11}
KRAKEN_REPORT_SUFFIX=${12}
HUMAN_ALIGNMENT_METRICS_SUFFIX=${13}
CONTIG_ALIGNMENT_METRICS_SUFFIX=${14}
NCBI_BLAST_TOOL_DIR=${15}
MASK_LOW_COMPLEXITY=${16}
TOOLS_DIR=${17}

echo -e "START: $(date)\nBacteria Pipeline\nFastq dir: $FASTQ_DIR\nResults dir: $RESULTS_DIR"
cd $RESULTS_DIR

ml python/3.6.1 java 
ml biology bwa samtools gatk

SAMPLE=$(find ${FASTQ_DIR}/ -maxdepth 1 -name "*${R1_SUFFIX}" -exec basename {} \; | \
    grep -v "Undetermined" | sed "s/${R1_SUFFIX}//" | sed -n ${SLURM_ARRAY_TASK_ID}p)
export PATH=${KRAKEN_TOOL_DIR}:$PATH

echo "Sample: $SAMPLE"

echo "### Trimming fastqs ### - START: $(date)"
java -jar ${TOOLS_DIR}/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 -trimlog \
    ${SAMPLE}_trimmomatic_log.txt ${FASTQ_DIR}/${SAMPLE}${R1_SUFFIX} ${FASTQ_DIR}/${SAMPLE}${R2_SUFFIX} \
    ${SAMPLE}_trimmed${R1_SUFFIX} ${SAMPLE}_unpaired_trimmed${R1_SUFFIX} \
    ${SAMPLE}_trimmed${R2_SUFFIX} ${SAMPLE}_unpaired_trimmed${R2_SUFFIX} \
    ILLUMINACLIP:${TOOLS_DIR}/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 MINLEN:36
echo "### Trimming fastqs ### - END: $(date)"

echo "### Counting fastq read counts ### - START: $(date)"
PRE_TRIM_READ_COUNT=$(echo $(zcat ${FASTQ_DIR}/${SAMPLE}${R1_SUFFIX} | wc -l ) \
    $(zcat ${FASTQ_DIR}/${SAMPLE}${R2_SUFFIX} | wc -l) | awk '{ print ($1 + $2) / 4 }' )
PAIR_TRIM_READ_COUNT=$(echo $(zcat ${SAMPLE}_trimmed${R1_SUFFIX} | wc -l ) \
    $(zcat ${SAMPLE}_trimmed${R2_SUFFIX} | wc -l) | awk '{ print ($1 + $2) / 4 }' )
UNPAIR_TRIM_READ_COUNT=$(echo $(zcat ${SAMPLE}_unpaired_trimmed${R1_SUFFIX} | wc -l ) \
    $(zcat ${SAMPLE}_unpaired_trimmed${R2_SUFFIX} | wc -l) | awk '{ print ($1 + $2) / 4 }' )
echo -e "sample\tread_count\ttrimmed_read_count\tunpaired_trimmed_read_count" > ${SAMPLE}.read_counts.tsv
echo -e "$SAMPLE\t$PRE_TRIM_READ_COUNT\t$PAIR_TRIM_READ_COUNT\t$UNPAIR_TRIM_READ_COUNT" >> ${SAMPLE}.read_counts.tsv
echo "### Counting fastq read counts ### - END: $(date)"

echo "### Aligning sample to human ### - START: $(date)"
bwa aln -t 4 $REF_FASTA ${SAMPLE}_trimmed${R1_SUFFIX} > ${SAMPLE}_R1.sai
bwa aln -t 4 $REF_FASTA ${SAMPLE}_trimmed${R2_SUFFIX} > ${SAMPLE}_R2.sai
bwa sampe -a 700 $REF_FASTA ${SAMPLE}_R1.sai ${SAMPLE}_R2.sai \
    ${SAMPLE}_trimmed${R1_SUFFIX} ${SAMPLE}_trimmed${R2_SUFFIX} | \
    samtools view -b - | samtools sort -o ${SAMPLE}_human_aligned.bam -
samtools index ${SAMPLE}_human_aligned.bam
echo "### Aligning sample to human ### - START: $(date)"

echo "### Collecting human alignment metrics ### - START: $(date)"
gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" CollectAlignmentSummaryMetrics \
    -R $REF_FASTA -I ${SAMPLE}_human_aligned.bam -O ${SAMPLE}${HUMAN_ALIGNMENT_METRICS_SUFFIX}
echo "### Collecting human alignment metrics ### - END: $(date)"

echo "### Getting meta data for human matches from BAM ### - START: $(date)"
samtools view  ${SAMPLE}_human_aligned.bam $(samtools view -H ${SAMPLE}_human_aligned.bam | grep -P "@SQ.+SN:chr" | \
    cut -f 2 | sed 's/^SN://')| cut -f 1 | sort -u > ${SAMPLE}_any_mapping_to_human_query_names.txt
# samtools view ${SAMPLE}_human_aligned.bam $(samtools view -H ${SAMPLE}_human_aligned.bam | grep -P "@SQ.+SN:chr..?\t" | \
#     cut -f 2 | sed 's/^SN://') | awk '$7=="="' | cut -f 1 | sort -u > ${SAMPLE}_both_mapping_to_human_query_names.txt
echo "### Getting meta data for human matches from BAM ### - END: $(date)"

echo "### Filtering BAM for reads that don't match human ### - START: $(date)"
gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" FilterSamReads \
    --FILTER excludeReadList -I ${SAMPLE}_human_aligned.bam -O ${SAMPLE}_no_human.bam \
    -RLF ${SAMPLE}_any_mapping_to_human_query_names.txt --VALIDATION_STRINGENCY SILENT
gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" SamToFastq -I ${SAMPLE}_no_human.bam \
    -F ${SAMPLE}_no_human${R1_SUFFIX} -F2 ${SAMPLE}_no_human${R2_SUFFIX} --VALIDATION_STRINGENCY SILENT
echo "### Filtering BAM for reads that don't match human ### - END: $(date)"

echo "### De novo assembling contigs ### - START: $(date)"
mkdir contigs_${SAMPLE}
python3 /oak/stanford/groups/cgawad/Sequencing_Analysis_Tools/SPAdes-3.14.0-Linux/bin/spades.py \
    -t 4 -m 64 -1 ${SAMPLE}_no_human${R1_SUFFIX} -2 ${SAMPLE}_no_human${R2_SUFFIX} -o contigs_${SAMPLE}
mv contigs_${SAMPLE}/contigs.fasta ${SAMPLE}${CONTIGS_SUFFIX}
echo "### De novo assembling contigs ### - END: $(date)"

echo "### Aligning reads to contigs ### - START: $(date)"
bwa index ${SAMPLE}${CONTIGS_SUFFIX}
bwa mem -M ${SAMPLE}${CONTIGS_SUFFIX} ${SAMPLE}_no_human${R1_SUFFIX} ${SAMPLE}_no_human${R2_SUFFIX} | \
    samtools view -b - | samtools sort -o ${SAMPLE}_contig_aligned.bam -
echo "### Aligning reads to contigs ### - END: $(date)"

echo "### Collecting contig alignment metrics ### - START: $(date)"
gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" CollectAlignmentSummaryMetrics \
    -R ${SAMPLE}${CONTIGS_SUFFIX} -I ${SAMPLE}_contig_aligned.bam -O ${SAMPLE}${CONTIG_ALIGNMENT_METRICS_SUFFIX}
echo "### Collecting contig alignment metrics ### - END: $(date)"

echo "### Exporting summarized and contig read targets from BAM ### - START: $(date)"
HUMAN_ALIGNED_READS=$(samtools view ${SAMPLE}_human_aligned.bam | cut -f 3 | grep "chr" | wc -l)
CONTIG_ALIGNED_READS=$(samtools view ${SAMPLE}_contig_aligned.bam | cut -f 3 | grep "NODE" | wc -l)
UNALIGNED_READS=$(samtools view ${SAMPLE}_contig_aligned.bam | cut -f 3 | grep -v "NODE" | wc -l)
echo -e "sample\thuman_aligned\tcontig_aligned\tunaligned" > ${SAMPLE}${SUMMED_READ_TARGETS_SUFFIX}
echo -e "${SAMPLE}\t${HUMAN_ALIGNED_READS}\t${CONTIG_ALIGNED_READS}\t${UNALIGNED_READS}" >> ${SAMPLE}${SUMMED_READ_TARGETS_SUFFIX}

echo -e "sample\ttarget" > ${SAMPLE}${CONTIGS_READ_TARGETS_SUFFIX}
samtools view ${SAMPLE}_contig_aligned.bam | cut -f 3 | grep "NODE" > ${SAMPLE}_temp_contig_read_targets.txt
printf "${SAMPLE}\n%0.s" $(seq $(cat ${SAMPLE}_temp_contig_read_targets.txt | wc -l)) | \
    paste - ${SAMPLE}_temp_contig_read_targets.txt >> ${SAMPLE}${CONTIGS_READ_TARGETS_SUFFIX}
echo "### Exporting summarized and contig read targets from BAM ### - END: $(date)"

echo "### Running kraken2 on non-human matches ### - START: $(date)"
for DB_TYPE in ${KRAKEN_DB_TYPES_ARRAY[@]}; do
    kraken2 --db ${KRAKEN_DB_DIR_PREFIX}${DB_TYPE} --threads 4 --output ${SAMPLE}_no_human_vs_kraken.tsv --paired --gzip-compressed \
        --report ${SAMPLE}_${DB_TYPE}${KRAKEN_REPORT_SUFFIX} ${SAMPLE}_no_human${R1_SUFFIX} ${SAMPLE}_no_human${R2_SUFFIX}
done
echo "### Running kraken2 on non-human matches ### - END: $(date)"

if [ $MASK_LOW_COMPLEXITY -eq 1 ]; then
    echo "### Marking low complexity regions in contigs ### - START: $(date)"
    mv ${SAMPLE}${CONTIGS_SUFFIX} ${SAMPLE}${CONTIGS_SUFFIX}_temp
    ${NCBI_BLAST_TOOL_DIR}/bin/dustmasker -in ${SAMPLE}${CONTIGS_SUFFIX}_temp -outfmt fasta -out ${SAMPLE}${CONTIGS_SUFFIX}
    rm ${SAMPLE}${CONTIGS_SUFFIX}_temp 
    echo "### Marking low complexity regions in contigs ### - END: $(date)"
fi

if [ ! -f ${SAMPLE}_no_human_vs_kraken.tsv ]; then
    echo "Final file ${SAMPLE}_no_human_vs_kraken.tsv not found. Exiting with code 1"
    exit 1
fi
rm ${SAMPLE}_trimmomatic_log.txt
rm ${SAMPLE}_trimmed${R1_SUFFIX} ${SAMPLE}_trimmed${R2_SUFFIX}
rm ${SAMPLE}_unpaired_trimmed${R1_SUFFIX} ${SAMPLE}_unpaired_trimmed${R2_SUFFIX}
rm -r contigs_${SAMPLE} 
rm ${SAMPLE}${CONTIGS_SUFFIX}.amb ${SAMPLE}${CONTIGS_SUFFIX}.ann ${SAMPLE}${CONTIGS_SUFFIX}.bwt
rm ${SAMPLE}${CONTIGS_SUFFIX}.pac ${SAMPLE}${CONTIGS_SUFFIX}.sa
rm ${SAMPLE}_R1.sai ${SAMPLE}_R2.sai
rm ${SAMPLE}_human_aligned.bam* ${SAMPLE}_no_human.bam* ${SAMPLE}_contig_aligned.bam*
rm ${SAMPLE}_temp_contig_read_targets.txt
rm ${SAMPLE}_no_human${R1_SUFFIX} ${SAMPLE}_no_human${R2_SUFFIX}
rm ${SAMPLE}_any_mapping_to_human_query_names.txt
rm ${SAMPLE}_no_human_vs_kraken.tsv
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"