import pandas as pd
import numpy as np
import datetime
import argparse
import glob
import json
import sys
import re


parser = argparse.ArgumentParser(
    description="""Converts a BLAST JSON file into a TSV""")
parser.add_argument('-p', '--project', help="Current project name", required=True)
parser.add_argument('-k', '--kraken', help="Kraken DB type", default="microbial")
parser.add_argument('-d', '--identify', help="Current identify name", default="")
parser.add_argument('-o', '--suffix', help="Suffix for output", default="")
parser.add_argument('-i', '--input', help="Input kraken reports file", default="")
parser.add_argument('-s', '--samples', help="String of sample names", default="")
parser.add_argument('-q', '--query', help="Query", default="contig")
parser.add_argument('-g', '--add_genus', help="Add genus", action="store_true", default=False)
args = parser.parse_args()


if len(args.input) == 0:
    args.input = "{}.{}.kraken_reports.tsv".format(args.project, args.kraken)
if len(args.identify) == 0:
    args.identify = args.project
if len(args.suffix) == 0:
    args.suffix = ".{}_{}.kraken_jtree.json".format(args.kraken, args.query)
    if args.add_genus:
        args.suffix = ".genus_{}.kraken_jtree.json".format(args.query)
if len(args.samples) > 0:
    samples = args.samples.split(":")
if len(args.samples) == 0:
    samples = glob.glob("*_{}s.fasta".format(args.query))
    for i in range(len(samples)):
        samples[i] = samples[i].replace("*_{}s.fasta".format(args.query), '')


kraken_df = pd.read_csv(args.input, sep="\t")
# samples = kraken_df['sample'].unique().tolist()
for sample in samples:
    if sample.count("_") == 6:
        short_sample = sample.split("_")[5]
    else:
        short_sample = sample
    jtree_data = []
    sample_df = kraken_df.loc[kraken_df['sample'] == sample]
    current_newick = ''
    saved_newicks = pd.DataFrame()
    previous_spaces = 0
    increment = 0
    unique_depths = []
    for row in reversed(sample_df.index):
        if len(sample_df['sciname'][row].strip().split()) > 2:
            continue
        if args.add_genus and len(sample_df['sciname'][row].strip().split()) > 1:
            continue
        increment += 1
        new_spaces = sample_df['sciname'][row].count(' ')
        new_newick = '%s:%s{%s}' % (re.sub(r'[^A-Za-z0-9 ]+', '', sample_df['sciname'][row].strip()), str(new_spaces), str(increment))
        new_dict = {
            'edge_num': int(increment),
            'percent_fragments_covered': float(sample_df['percent_fragments_covered'][row]),
            'fragments_covered': int(sample_df['fragments_covered'][row]),
            'fragments_assigned': int(sample_df['fragments_assigned'][row]),
            'rank_code': str(sample_df['rank_code'][row]),
            'taxid': str(sample_df['taxid'][row])
        }
        jtree_data.append(new_dict)
        if new_spaces not in unique_depths:
            unique_depths.append(new_spaces)
        if len(saved_newicks) > 0:
            if len(saved_newicks.loc[saved_newicks.iloc[:,1] == new_spaces]) > 0:
                previous_newicks = ",".join(saved_newicks.loc[saved_newicks.iloc[:, 1] == new_spaces][0])
                new_newick = '{},{}'.format(new_newick, previous_newicks)
                saved_newicks = saved_newicks[saved_newicks.iloc[:, 1] != new_spaces]
        if len(current_newick) < 1:
            current_newick = new_newick
        elif new_spaces < previous_spaces:
            current_newick = '({}){}'.format(current_newick, new_newick)
        elif new_spaces == previous_spaces:
            current_newick = '{},{}'.format(current_newick, new_newick)
        elif new_spaces > previous_spaces:
            saved_newicks = saved_newicks.append(pd.Series([current_newick, previous_spaces]), ignore_index=True)
            current_newick = new_newick
        previous_spaces = new_spaces
    if "unclassified" not in current_newick:
        current_newick = '({},);'.format(current_newick)
    else:
        current_newick = '({});'.format(current_newick)
    jtree_dict = {
        'tree': str(current_newick),
        'data': jtree_data,
        'metadata': {'info': 'kraken2 report to json', 
                     'data': int(len(unique_depths)), 
                     'max_depth': int(len(unique_depths))}
    }
    with open('{}.{}{}'.format(args.project, short_sample, args.suffix), 'w') as outfile:
        json.dump(jtree_dict, outfile)