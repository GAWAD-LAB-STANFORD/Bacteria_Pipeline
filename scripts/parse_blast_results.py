import pandas as pd
import json
import glob
import sys
import re

db_type = sys.argv[1]
blast_results_prefix = sys.argv[2]
parsed_blast_results_filename = sys.argv[3]
blastn_jsons = glob.glob("{}*.json".format(blast_results_prefix))
print("{} blast json result files found with {} prefix".format(len(blastn_jsons), blast_results_prefix))
all_top_blast_hits_list_list = []

# for each blast result file, consolidate the top results
for blastn_json in blastn_jsons:
    json_file = open(blastn_json, "r")
    blastn_dict = json.load(json_file)
    json_file.close()
    for blast_result in blastn_dict['BlastOutput2']:
        extended_contig_name = re.sub(" .+", "", blast_result['report']['results']['search']['query_title'])
        sample = extended_contig_name.split(".")[0]
        query_contig_name = ".".join(extended_contig_name.split(".")[1:])
        query_contig_length = blast_result['report']['results']['search']['query_len']
        hits = blast_result['report']['results']['search']['hits']
        if len(hits) > 0:
            for hit_count in range(len(hits)):
                hit =  hits[hit_count]
                if db_type == "plasmid":
                    before, keyword, after = hit['description'][0]['title'].partition("plasmid ")
                    plasmid = after.split(" ")[0][:-1]
                    source = " ".join(before.split(" ")[1:-1])
                    new_addition = [sample, hit_count + 1, query_contig_name, query_contig_length, source, plasmid, 
                                    hit['description'][0]['accession'], hit['hsps'][0]['align_len'], hit['len']]
                elif db_type == "viral":
                    before, keyword, after = hit['description'][0]['title'].partition("phage ")
                    phage = after.split(" ")[0][:-1]
                    source = " ".join(before.split(" ")[1:-1])
                    new_addition = [sample, hit_count + 1, query_contig_name, query_contig_length, source, phage, 
                                    hit['description'][0]['accession'], hit['hsps'][0]['align_len'], hit['len']]
                else:
                    new_addition = [sample, hit_count + 1, query_contig_name, query_contig_length, 
                                    hit['description'][0]['sciname'], hit['description'][0]['taxid'],
                                    hit['description'][0]['accession'], hit['hsps'][0]['align_len'], hit['len']]
                if '' not in new_addition:
                    all_top_blast_hits_list_list.append(new_addition)

if len(all_top_blast_hits_list_list) > 0:
    if db_type == "plasmid":
        blast_results_df = pd.DataFrame(all_top_blast_hits_list_list, 
                                    columns=['sample', 'hit_rank', 'contig', 'contig_length', 
                                    'source', 'plasmid', 'accession', 'top_hsp_align_len', 'reference_len'])
    elif db_type == "viral":
        blast_results_df = pd.DataFrame(all_top_blast_hits_list_list, 
                                    columns=['sample', 'hit_rank', 'contig', 'contig_length', 
                                    'source', 'phage', 'accession', 'top_hsp_align_len', 'reference_len'])
    else:
        blast_results_df = pd.DataFrame(all_top_blast_hits_list_list, 
                                        columns=['sample', 'hit_rank', 'contig', 'contig_length', 
                                        'species', 'hit_taxid', 'accession', 'top_hsp_align_len', 'reference_len'])     
    blast_results_df.to_csv(parsed_blast_results_filename, header = True, index = False, sep="\t")