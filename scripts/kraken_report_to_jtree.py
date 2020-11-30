import pandas as pd
import numpy as np
import datetime
import json
import sys

input_file_name = sys.argv[1]
project = sys.argv[2]
output_suffix = sys.argv[3]
samples = sys.argv[4].split(":")

kraken_df = pd.read_csv(input_file_name, sep="\t")
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
        increment += 1
        new_spaces = sample_df['sciname'][row].count(' ')
        new_newick = '%s:%s{%s}' % (sample_df['sciname'][row].strip(), str(new_spaces), str(increment))
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
    with open('{}.{}{}'.format(project, short_sample, output_suffix), 'w') as outfile:
        json.dump(jtree_dict, outfile)