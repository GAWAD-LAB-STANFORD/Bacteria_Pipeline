from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import glob
import sys
import os
import re

contig_file_suffix = sys.argv[1]
record_length_minimum = int(sys.argv[2])
records_per_file = int(sys.argv[3])
organized_contigs_prefix = sys.argv[4]
print("Finding contigs of at least {} bases long".format(record_length_minimum))

contigs = glob.glob("*{}".format(contig_file_suffix))
long_records = []
for contig in contigs:
    sample = re.sub(contig_file_suffix, "", contig)
    long_contig_count = 0
    contigs_with_low_complexity = 0
    sample_records = list(SeqIO.parse(contig, "fasta"))
    for record in sample_records:
        if len(record) >= record_length_minimum:
            long_contig_count += 1
            record.id = "{}.{}".format(sample, record.id)
            long_records.append(record)
    print("\t{} long contigs found for sample {}".format(long_contig_count, sample))
    print("\t{} contigs had low complexity regions that were removed".format(contigs_with_low_complexity))
print("{} total long contigs found across all samples".format(len(long_records)))

file_count = 0
if len(long_records) <= records_per_file:
    file_count += 1
    SeqIO.write(long_records, "{}{}.fasta".format(organized_contigs_prefix, 1), "fasta")
else:
    for end_index in range(records_per_file, len(long_records), records_per_file):
        file_count += 1
        SeqIO.write(long_records[end_index-records_per_file:end_index],
            "{}{:.0f}.fasta".format(organized_contigs_prefix, file_count), "fasta")
print("Split into {} files with {} records each".format(file_count, records_per_file))
        