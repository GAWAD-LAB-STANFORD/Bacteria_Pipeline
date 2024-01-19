import pandas as pd
import json
import glob
import sys
import re


db_type = sys.argv[1]
blast_results_prefix = sys.argv[2]
parsed_blast_results_filename = sys.argv[3]
pipeline_status_filename = sys.argv[4]
blastn_jsons = glob.glob("{}*.json".format(blast_results_prefix))
pipeline_status_output = "{} blast json result files found with {} prefix\n".format(len(blastn_jsons), blast_results_prefix)
all_top_blast_hits_list_list = []


# for each blast result file, consolidate the top results
for blastn_json in blastn_jsons:
    json_file = open(blastn_json, "r")
    blastn_dict = json.load(json_file)
    json_file.close()
    for blast_result in blastn_dict['BlastOutput2']:
        extended_query_name = re.sub(" .+", "", blast_result['report']['results']['search']['query_title'])
        sample = extended_query_name.split(".")[0]
        query_name = ".".join(extended_query_name.split(".")[1:])
        query_length = blast_result['report']['results']['search']['query_len']
        hits = blast_result['report']['results']['search']['hits']
        if len(hits) > 0:
            for hit_count in range(len(hits)):
                hit =  hits[hit_count]
                fraction_identical = round(hit['hsps'][0]['identity'] / query_length, 4)
                if db_type == "plasmid":
                    before, keyword, after = hit['description'][0]['title'].partition("plasmid ")
                    plasmid = after.split(" ")[0][:-1]
                    source = " ".join(before.split(" ")[1:-1])
                    new_addition = [sample, hit_count + 1, query_name, query_length,
                                    hit['description'][0]['title'], hit['description'][0]['accession'],
                                    plasmid, source, hit['hsps'][0]['align_len'], hit['len'],
                                    hit['hsps'][0]['query_from'], hit['hsps'][0]['query_to'],
                                    hit['hsps'][0]['hit_from'], hit['hsps'][0]['hit_to'],
                                    hit['hsps'][0]['gaps'], hit['hsps'][0]['evalue'], fraction_identical]
                elif db_type == "viral":
                    before, keyword, after = hit['description'][0]['title'].partition("phage ")
                    phage = after.split(" ")[0][:-1]
                    source = " ".join(before.split(" ")[1:-1])
                    new_addition = [sample, hit_count + 1, query_name, query_length, 
                                    hit['description'][0]['title'], hit['description'][0]['accession'],
                                    phage, source, hit['hsps'][0]['align_len'], hit['len'],
                                    hit['hsps'][0]['query_from'], hit['hsps'][0]['query_to'],
                                    hit['hsps'][0]['hit_from'], hit['hsps'][0]['hit_to'],
                                    hit['hsps'][0]['gaps'], fraction_identical]
                else:
                    new_addition = [sample, hit_count + 1, query_name, query_length, 
                                    hit['description'][0]['sciname'], hit['description'][0]['taxid'],
                                    hit['description'][0]['accession'], hit['hsps'][0]['align_len'], hit['len'],
                                    hit['hsps'][0]['query_from'], hit['hsps'][0]['query_to'],
                                    hit['hsps'][0]['hit_from'], hit['hsps'][0]['hit_to'],
                                    hit['hsps'][0]['gaps'], fraction_identical]
                all_top_blast_hits_list_list.append(new_addition)


if len(all_top_blast_hits_list_list) > 0:
    if db_type == "plasmid":
        blast_results_df = pd.DataFrame(all_top_blast_hits_list_list, 
                                        columns=['sample', 'hit_rank', 'query', 'query_length', 'full_name',
                                        'accession', 'plasmid', 'source', 'top_hsp_align_len', 'reference_len',
                                        'query_from', 'query_to', 'hit_from', 'hit_to', 'gaps', 'evalue', 'fraction_identical'])
    elif db_type == "viral":
        blast_results_df = pd.DataFrame(all_top_blast_hits_list_list, 
                                        columns=['sample', 'hit_rank', 'query', 'query_length', 'full_name',
                                        'accession', 'virus', 'source', 'top_hsp_align_len', 'reference_len',
                                        'query_from', 'query_to', 'hit_from', 'hit_to', 'gaps', 'evalue', 'fraction_identical'])
    else:
        blast_results_df = pd.DataFrame(all_top_blast_hits_list_list, 
                                        columns=['sample', 'hit_rank', 'query', 'query_length', 
                                        'species', 'hit_taxid', 'accession', 'top_hsp_align_len', 'reference_len',
                                        'query_from', 'query_to', 'hit_from', 'hit_to', 'gaps', 'evalue', 'fraction_identical'])     
    blast_results_df.to_csv(parsed_blast_results_filename, header = True, index = False, sep="\t")


file1 = open(pipeline_status_filename, "a")
file1.write(pipeline_status_output) 
file1.close()