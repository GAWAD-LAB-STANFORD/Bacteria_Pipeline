#!/bin/bash
PIPELINE_STATUS=$1
PROJECT=$2
KRAKEN_DB_TYPES=$3
RUN_DIR=$4
SAMPLE_SHEET=$5

echo "### Summarizing metrics ### - START: $(date)" >> $PIPELINE_STATUS
SAMPLE_READ_COUNTS="${PROJECT}.sample_read_counts.tsv"
READ_COUNT_FILENAMES=( $(ls *_read_counts.tsv) )
head -n 1 ${READ_COUNT_FILENAMES[0]} > $SAMPLE_READ_COUNTS
for i in ${READ_COUNT_FILENAMES[@]}; do tail -n +2 $i; done >> $SAMPLE_READ_COUNTS
if [ ! -z $RUN_DIR ]; then
    DESIRED_CLUSTER_DENSITY=$(sed -n $(grep -n "Lane" $SAMPLE_SHEET | sed "s/:Lane.*//" | \
        awk '{print $1 + 1}')p $SAMPLE_SHEET | cut -d , -f $(grep "Lane" $SAMPLE_SHEET | \
        sed "s/,/\n/g" | nl | grep "Desired_Cluster_Density" | cut -f 1))
    INITAL_CONCENTRATION_COLUMN=$(grep "Initial_Concentration" $SAMPLE_SHEET | tr , "\n" | nl | grep "Initial_Concentration" | cut -f 1)
    LIBRARY_GROUP_COLUMN=$(grep "Library_Group" $SAMPLE_SHEET | tr , "\n" | nl | grep "Library_Group" | cut -f 1)
    if [ ! -z $INITAL_CONCENTRATION_COLUMN ] && [ ! -z $DESIRED_CLUSTER_DENSITY ]; then
        echo "### Calculating library concentration corrections ### - START: $(date)" >> $PIPELINE_STATUS
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
        echo "### Calculating library concentration corrections ### - END: $(date)" >> $PIPELINE_STATUS
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
echo "Merged human alignment metrics" >> $PIPELINE_STATUS

CONTIG_ALIGNMENT_METRICS_FILENAMES=( $(ls *_contig_alignment_metrics.tsv) )
echo -e sample"\t"$(head -n 7 ${CONTIG_ALIGNMENT_METRICS_FILENAMES[0]} | tail -n 1) | sed 's/ /\t/g' > ${PROJECT}.contig_alignment_metrics.tsv
for i in ${CONTIG_ALIGNMENT_METRICS_FILENAMES[@]}; do 
    SAMPLE=$(echo $i | sed "s/_contig_alignment_metrics.tsv//")
    R1=$(head -n 8 $i | tail -n 1)
    R2=$(head -n 9 $i | tail -n 1)
    PAIR=$(head -n 10 $i | tail -n 1)
    echo -e "$SAMPLE\t$R1\n$SAMPLE\t$R2\n$SAMPLE\t$PAIR"
done | sed 's/ /\t/g' >> ${PROJECT}.contig_alignment_metrics.tsv
echo "Merged contig alignment metrics" >> $PIPELINE_STATUS

KRAKEN_DB_TYPES_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
for DB_TYPE in ${KRAKEN_DB_TYPES_ARRAY[@]}; do
    echo -e "sample\tpercent_fragments_covered\tfragments_covered\tfragments_assigned\trank_code\ttaxid\tsciname" > \
        ${PROJECT}.${DB_TYPE}.kraken_reports.tsv
    KRAKEN_REPORT_FILENAMES=( $(ls *_${DB_TYPE}_kraken_report.tsv) )
    for i in ${KRAKEN_REPORT_FILENAMES[@]}; do
        SAMPLE=$(echo $i | sed "s/_${DB_TYPE}_kraken_report.tsv//")
        cat $i | sed 's/^ \+/'${SAMPLE}'\t/' | awk '$2>=1.00' >> ${PROJECT}.${DB_TYPE}.kraken_reports.tsv
    done
done
KRAKEN_REPORT_FILENAMES=( $(ls *_kraken_report.tsv) )
echo "Merged kraken reports" >> $PIPELINE_STATUS

echo -e "sample\tfasta_header" > ${PROJECT}.contig_data.tsv
CONTIGS_FILENAMES=( $(ls *_contigs.fasta) )
for i in ${CONTIGS_FILENAMES[@]}; do 
    grep ">" $i | xargs -i echo -e $(echo $i | sed "s/_contigs.fasta//")"\t"{} >> ${PROJECT}.contig_data.tsv
done
echo "Merged contigs" >> $PIPELINE_STATUS

SUMMED_READ_TARGETS_FILENAMES=( $(ls *_summed_read_targets.tsv) )
head -n 1 ${SUMMED_READ_TARGETS_FILENAMES[0]} > ${PROJECT}.summed_read_targets.tsv
for i in ${SUMMED_READ_TARGETS_FILENAMES[@]}; do sed -n 2p $i >> ${PROJECT}.summed_read_targets.tsv; done
echo "Merged summarized read targets" >> $PIPELINE_STATUS

CONTIG_READ_TARGETS_FILENAMES=( $(ls *_contig_read_targets.tsv) )
head -n 1 ${CONTIG_READ_TARGETS_FILENAMES[0]} > ${PROJECT}.contig_read_targets.tsv
for i in ${CONTIG_READ_TARGETS_FILENAMES[@]}; do tail -n +2 $i >> ${PROJECT}.contig_read_targets.tsv; done
echo "Merged contig read targets" >> $PIPELINE_STATUS

rm ${READ_COUNT_FILENAMES[@]} ${KRAKEN_REPORT_FILENAMES[@]}
rm ${SUMMED_READ_TARGETS_FILENAMES[@]} ${CONTIG_READ_TARGETS_FILENAMES[@]}
rm ${HUMAN_ALIGNMENT_METRICS_FILENAMES[@]} ${CONTIG_ALIGNMENT_METRICS_FILENAMES[@]}
echo "Deleted intermediate files" >> $PIPELINE_STATUS
echo "### Summarizing metrics ### - END: $(date)" >> $PIPELINE_STATUS