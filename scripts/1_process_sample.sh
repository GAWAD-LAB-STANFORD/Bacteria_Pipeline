#!/bin/bash
#
#SBATCH --job-name=1_process_sample
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
FASTQ_DIR=$1
RESULTS_DIR=$2
R1_SUFFIX=$3
R2_SUFFIX=$4
SKIP_TRIMMOMATIC=$5
REF_FASTA_ARRAY=( $(echo $6 | sed 's/:/ /g') )
REF_NAME_ARRAY=( $(echo $7 | sed 's/:/ /g') )
KRAKEN_DB_TYPE_ARRAY=( $(echo $8 | sed 's/-/ /g') )
KRAKEN_DB_DIR_PREFIX=$9
TOOLS_DIR=${10}
SAMPLE_ARRAY=( $(echo ${11} | sed 's/:/ /g') )
SAMPLE=${SAMPLE_ARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}

echo -e "START: $(date)\nBacteria Pipeline\nFastq dir: $FASTQ_DIR\nResults dir: $RESULTS_DIR\nSample: $SAMPLE"
cd $RESULTS_DIR

ml python/3.6.1 java 
ml biology bwa samtools gatk

R1_FASTQ=${FASTQ_DIR}/${SAMPLE}${R1_SUFFIX}
R2_FASTQ=${FASTQ_DIR}/${SAMPLE}${R2_SUFFIX}
export PATH=${TOOLS_DIR}/kraken2-2.0.8-beta:$PATH

if [ $SKIP_TRIMMOMATIC -eq 0 ]; then
    UNTRIMMED_R1_FASTQ=$R1_FASTQ
    UNTRIMMED_R2_FASTQ=$R2_FASTQ
    R1_FASTQ=$(echo ${SAMPLE}${R1_SUFFIX} | sed "s/_R1_/_R1_trimmed_/")
    R2_FASTQ=$(echo ${SAMPLE}${R2_SUFFIX} | sed "s/_R2_/_R2_trimmed_/")
    UNPAIRED_R1_FASTQ=$(echo ${SAMPLE}${R1_SUFFIX} | sed "s/_R1_/_R1_trimmed_unpaired_/")
    UNPAIRED_R2_FASTQ=$(echo ${SAMPLE}${R2_SUFFIX} | sed "s/_R2_/_R2_trimmed_unpaired_/")
    
    echo "### Trimming fastqs ### - START: $(date)"
    java -jar ${TOOLS_DIR}/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 -trimlog \
        ${SAMPLE}_trimmomatic_log.txt ${UNTRIMMED_R1_FASTQ} ${UNTRIMMED_R2_FASTQ} \
        ${R1_FASTQ} ${UNPAIRED_R1_FASTQ} \
        ${R2_FASTQ} ${UNPAIRED_R2_FASTQ} \
        ILLUMINACLIP:${TOOLS_DIR}/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 MINLEN:36
    echo "### Trimming fastqs ### - END: $(date)"
fi

echo "### Counting fastq read counts ### - START: $(date)"
READ_COUNT=$(echo $(zcat $R1_FASTQ | wc -l ) \
    $(zcat $R2_FASTQ | wc -l) | awk '{ print ($1 + $2) / 4 }' )
echo -e "sample\tread_count" > ${SAMPLE}.read_counts.tsv
echo -e "$SAMPLE\t$READ_COUNT" >> ${SAMPLE}.read_counts.tsv
echo "### Counting fastq read counts ### - END: $(date)"

PREV_R1_FASTQ=$R1_FASTQ
PREV_R2_FASTQ=$R2_FASTQ
for ((REF_INDEX = 0 ; REF_INDEX < ${#REF_FASTA_ARRAY[@]} ; REF_INDEX++)); do
    REF_FASTA=${REF_FASTA_ARRAY[$REF_INDEX]}
    REF_NAME=${REF_NAME_ARRAY[$REF_INDEX]}
    
    echo "### Aligning sample to $REF_NAME ### - START: $(date)"
    echo -e "Ref fasta: $REF_FASTA\nR1 fastq: $PREV_R1_FASTQ\nR2 fastq: $PREV_R2_FASTQ"
    bwa aln -t 4 $REF_FASTA $PREV_R1_FASTQ > ${SAMPLE}_${REF_NAME}_R1.sai
    bwa aln -t 4 $REF_FASTA $PREV_R2_FASTQ > ${SAMPLE}_${REF_NAME}_R2.sai
    bwa sampe -a 700 $REF_FASTA ${SAMPLE}_${REF_NAME}_R1.sai ${SAMPLE}_${REF_NAME}_R2.sai \
        $PREV_R1_FASTQ $PREV_R2_FASTQ | samtools view -b - | samtools sort -o ${SAMPLE}_${REF_NAME}_aligned.bam -
    samtools index ${SAMPLE}_${REF_NAME}_aligned.bam
    echo "### Aligning sample to $REF_NAME ### - END: $(date)"
    
    echo "### Collecting $REF_NAME alignment metrics ### - START: $(date)"
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" CollectAlignmentSummaryMetrics \
        -R $REF_FASTA -I ${SAMPLE}_${REF_NAME}_aligned.bam -O ${SAMPLE}_${REF_NAME}_alignment_metrics.tsv
    echo "### Collecting $REF_NAME alignment metrics ### - END: $(date)"
    
    echo "### Getting meta data for $REF_NAME matches from BAM ### - START: $(date)"
    samtools view  ${SAMPLE}_${REF_NAME}_aligned.bam $(samtools view -H ${SAMPLE}_${REF_NAME}_aligned.bam | grep -P "@SQ.+SN:chr" | \
        cut -f 2 | sed 's/^SN://')| cut -f 1 | sort -u > ${SAMPLE}_any_mapping_to_${REF_NAME}_query_names.txt
    # samtools view ${SAMPLE}_${REF_NAME}_aligned.bam $(samtools view -H ${SAMPLE}_${REF_NAME}_aligned.bam | grep -P "@SQ.+SN:chr..?\t" | \
    #     cut -f 2 | sed 's/^SN://') | awk '$7=="="' | cut -f 1 | sort -u > ${SAMPLE}_both_mapping_to_${REF_NAME}_query_names.txt
    echo "### Getting meta data for $REF_NAME matches from BAM ### - END: $(date)"

    echo "### Filtering BAM for reads that don't match $REF_NAME ### - START: $(date)"
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" FilterSamReads \
        --FILTER excludeReadList -I ${SAMPLE}_${REF_NAME}_aligned.bam -O ${SAMPLE}_no_${REF_NAME}.bam \
        -RLF ${SAMPLE}_any_mapping_to_${REF_NAME}_query_names.txt --VALIDATION_STRINGENCY SILENT
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" SamToFastq -I ${SAMPLE}_no_${REF_NAME}.bam \
        -F ${SAMPLE}_no_${REF_NAME}${R1_SUFFIX} -F2 ${SAMPLE}_no_${REF_NAME}${R2_SUFFIX} --VALIDATION_STRINGENCY SILENT
    echo "### Filtering BAM for reads that don't match $REF_NAME ### - END: $(date)"
    
    PREV_R1_FASTQ=${SAMPLE}_no_${REF_NAME}${R1_SUFFIX}
    PREV_R2_FASTQ=${SAMPLE}_no_${REF_NAME}${R2_SUFFIX}
    rm ${SAMPLE}_${REF_NAME}_R1.sai ${SAMPLE}_${REF_NAME}_R2.sai
    rm ${SAMPLE}_any_mapping_to_${REF_NAME}_query_names.txt
    # rm ${SAMPLE}_no_${REF_NAME}.bam ${SAMPLE}_no_${REF_NAME}.bam.bai
done

for ((REF_INDEX = 0 ; REF_INDEX < ${#REF_FASTA_ARRAY[@]} ; REF_INDEX++)); do
    REF_FASTA=${REF_FASTA_ARRAY[$REF_INDEX]}
    REF_NAME=${REF_NAME_ARRAY[$REF_INDEX]}
    
    NEXT_INDEX=$((REF_INDEX+1))
    if [ $NEXT_INDEX -eq ${#REF_FASTA_ARRAY[@]} ]; then
        mv ${SAMPLE}_no_${REF_NAME}${R1_SUFFIX} ${SAMPLE}_ref_filtered${R1_SUFFIX}
        mv ${SAMPLE}_no_${REF_NAME}${R2_SUFFIX} ${SAMPLE}_ref_filtered${R2_SUFFIX}
    else
        echo "not removing filtered fastqs"
        # rm ${SAMPLE}_no_${REF_NAME}${R1_SUFFIX} ${SAMPLE}_no_${REF_NAME}${R2_SUFFIX}
    fi
done

echo "### De novo assembling contigs ### - START: $(date)"
mkdir contigs_${SAMPLE}
python3 /oak/stanford/groups/cgawad/Sequencing_Analysis_Tools/SPAdes-3.14.0-Linux/bin/spades.py \
    -t 4 -m 64 -1 ${SAMPLE}_ref_filtered${R1_SUFFIX} -2 ${SAMPLE}_ref_filtered${R2_SUFFIX} -o contigs_${SAMPLE}
mv contigs_${SAMPLE}/contigs.fasta ${SAMPLE}_contigs.fasta
echo "### De novo assembling contigs ### - END: $(date)"

echo "### Aligning reads to contigs ### - START: $(date)"
bwa index ${SAMPLE}_contigs.fasta
bwa mem -M ${SAMPLE}_contigs.fasta ${SAMPLE}_ref_filtered${R1_SUFFIX} ${SAMPLE}_ref_filtered${R2_SUFFIX} | \
    samtools view -b - | samtools sort -o ${SAMPLE}_contig_aligned.bam -
echo "### Aligning reads to contigs ### - END: $(date)"

echo "### Collecting contig alignment metrics ### - START: $(date)"
gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" CollectAlignmentSummaryMetrics \
    -R ${SAMPLE}_contigs.fasta -I ${SAMPLE}_contig_aligned.bam -O ${SAMPLE}_contig_alignment_metrics.tsv
echo "### Collecting contig alignment metrics ### - END: $(date)"

echo "### Exporting summarized and contig read targets from BAM ### - START: $(date)"
READ_TARGETS_HEADER="sample"
READ_TARGETS_DATA="$SAMPLE"
for REF_NAME in ${REF_NAME_ARRAY[@]}; do
    TEMP_REF_ALIGNED_READS=$(samtools view ${SAMPLE}_${REF_NAME}_aligned.bam | cut -f 3 | grep "chr" | wc -l)
    READ_TARGETS_HEADER+="\t${REF_NAME}_aligned"
    READ_TARGETS_DATA+="\t${TEMP_REF_ALIGNED_READS}"
done
CONTIG_ALIGNED_READS=$(samtools view ${SAMPLE}_contig_aligned.bam | cut -f 3 | grep "NODE" | wc -l)
UNALIGNED_READS=$(samtools view ${SAMPLE}_contig_aligned.bam | cut -f 3 | grep -v "NODE" | wc -l)
echo -e "$READ_TARGETS_HEADER\tcontig_aligned\tunaligned" > ${SAMPLE}_summed_read_targets.tsv
echo -e "$READ_TARGETS_DATA\t$CONTIG_ALIGNED_READS\t$UNALIGNED_READS" >> ${SAMPLE}_summed_read_targets.tsv

echo -e "sample\ttarget" > ${SAMPLE}_contig_read_targets.tsv
samtools view ${SAMPLE}_contig_aligned.bam | cut -f 3 | grep "NODE" > ${SAMPLE}_temp_contig_read_targets.txt
printf "${SAMPLE}\n%0.s" $(seq $(cat ${SAMPLE}_temp_contig_read_targets.txt | wc -l)) | \
    paste - ${SAMPLE}_temp_contig_read_targets.txt >> ${SAMPLE}_contig_read_targets.tsv
echo "### Exporting summarized and contig read targets from BAM ### - END: $(date)"

# echo "### Marking low complexity regions in contigs ### - START: $(date)"
# ${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/dustmasker -in ${SAMPLE}_contigs.fasta -outfmt fasta -out ${SAMPLE}_contigs_high_complexity.fasta
# echo "### Marking low complexity regions in contigs ### - END: $(date)"

echo "### Running kraken2 on reference filtered matches ### - START: $(date)"
for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
    kraken2 --db ${KRAKEN_DB_DIR_PREFIX}${DB_TYPE} --threads 4 --output ${SAMPLE}_ref_filtered_vs_kraken.tsv --paired --gzip-compressed \
        --report ${SAMPLE}_${DB_TYPE}_kraken_report.tsv ${SAMPLE}_ref_filtered${R1_SUFFIX} ${SAMPLE}_ref_filtered${R2_SUFFIX}
done
echo "### Running kraken2 on reference filtered matches ### - END: $(date)"

if [ ! -f ${SAMPLE}_ref_filtered_vs_kraken.tsv ]; then
    echo "Final file ${SAMPLE}_ref_filtered_vs_kraken.tsv not found. Exiting with code 1"
    exit 1
fi
if [ $SKIP_TRIMMOMATIC -eq 0 ]; then
    rm ${SAMPLE}_trimmomatic_log.txt
    rm $R1_FASTQ $R2_FASTQ
    rm $UNPAIRED_R1_FASTQ $UNPAIRED_R2_FASTQ
fi
# for REF_NAME in ${REF_NAME_ARRAY[@]}; do
#     rm ${SAMPLE}_${REF_NAME}_aligned.bam ${SAMPLE}_${REF_NAME}_aligned.bam.bai
# done
rm -r contigs_${SAMPLE} 
rm ${SAMPLE}_contigs.fasta.amb ${SAMPLE}_contigs.fasta.ann ${SAMPLE}_contigs.fasta.bwt
rm ${SAMPLE}_contigs.fasta.pac ${SAMPLE}_contigs.fasta.sa
rm ${SAMPLE}_contig_aligned.bam*
rm ${SAMPLE}_temp_contig_read_targets.txt
# rm ${SAMPLE}_ref_filtered${R1_SUFFIX} ${SAMPLE}_ref_filtered${R2_SUFFIX}
rm ${SAMPLE}_any_mapping_to_ref_filtered_query_names.txt
rm ${SAMPLE}_ref_filtered_vs_kraken.tsv
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"