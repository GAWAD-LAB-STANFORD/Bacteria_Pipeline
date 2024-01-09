from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import glob
import sys
import os
import re


query_file_suffix = sys.argv[1]
record_length_minimum = int(sys.argv[2])
records_per_file = int(sys.argv[3])
organized_query_prefix = sys.argv[4]
pipeline_status_filename = sys.argv[5]


pipeline_status_output = "Finding queries of at least {} bases long\n".format(record_length_minimum)
queries = glob.glob("*{}".format(query_file_suffix))
long_records = []
for query in queries:
    sample = re.sub(query_file_suffix, "", query)
    long_query_count = 0
    sample_records = list(SeqIO.parse(query, "fasta"))
    for record in sample_records:
        if len(record) >= record_length_minimum:
            long_query_count += 1
            record.id = "{}.{}".format(sample, record.id)
            long_records.append(record)
    pipeline_status_output += "\t{} long queries found for sample {}\n".format(long_query_count, sample)


pipeline_status_output += "{} total long queriess found across all samples\n".format(len(long_records))
file_count = 0
if len(long_records) <= records_per_file:
    file_count += 1
    SeqIO.write(long_records, "{}{}.fasta".format(organized_query_prefix, 1), "fasta")
else:
    for end_index in range(records_per_file, len(long_records), records_per_file):
        file_count += 1
        SeqIO.write(long_records[end_index-records_per_file:end_index],
            "{}{:.0f}.fasta".format(organized_query_prefix, file_count), "fasta")
    if len(long_records) % records_per_file > 0:
        file_count += 1
        SeqIO.write(long_records[end_index:len(long_records)],
            "{}{:.0f}.fasta".format(organized_query_prefix, file_count), "fasta")


pipeline_status_output += "Split into {} files with {} records each\n".format(file_count, records_per_file)
file1 = open(pipeline_status_filename, "a")
file1.write(pipeline_status_output) 
file1.close()
        