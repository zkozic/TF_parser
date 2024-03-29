from Bio import SeqIO
import pandas as pd
import re
import sys
import warnings
warnings.filterwarnings("ignore")

#gb_file = sys.argv[1]
#gb_record = SeqIO.read(open(gb_file,"r"), "genbank")

def get_gene_id(desc_dict):
    if 'sym' in desc_dict:
        return desc_dict['sym']

    return desc_dict['acc']


gene_list = []
prom_starts_list = []
prom_ends_list = []
chr_list = []
strand_list = []
bs_starts_list = []
bs_ends_list = []
tf_list = []
scores_list = []
comp_list = []

for record in SeqIO.parse(open(gb_file,"r"), "genbank"):
    desc = record.description.split('|')
    desc_dict = {entry.split('=')[0].strip():entry.split('=')[1].strip() for entry in desc}

    gene = get_gene_id(desc_dict)
    chromosome = desc_dict['chr']
    strand = re.split("\(|\)", desc_dict['str'])[-2]
    prom_start = desc_dict['start']
    prom_end = desc_dict['end']

    for feature in record.features:
        feature = str(feature)
        fsplit = feature.split('\n')
        loc = re.split('\[|\]', [i for i in fsplit if 'location' in i][0])[1]
        tf = re.split("'|,", [i for i in fsplit if 'Value' in i][0])[2]
        bs_start = loc.split(':')[0]
        bs_end = loc.split(':')[1]
        score = re.split(" |'", [i for i in fsplit if 'Value' in i][0])[-2]
        comp = ''
        if '-' in feature:
            comp = 'Complement'
        gene_list += [gene]
        chr_list += [chromosome]
        strand_list += [strand]
        prom_starts_list += [prom_start]
        prom_ends_list += [prom_end]
        tf_list += [tf]
        bs_starts_list += [bs_start]
        bs_ends_list += [bs_end]
        scores_list += [score]
        comp_list += [comp]

result_df = list(zip(gene_list, chr_list, strand_list, prom_starts_list, prom_ends_list, tf_list, bs_starts_list, bs_ends_list, scores_list, comp_list))
result_df = pd.DataFrame(data = result_df, columns = ['Gene', 'Chr', 'Strand', 'Promoter start', 'Promoter end', 'Transcription factor', 'Binding site start', 'Binding site end', 'Score', 'Complement'])

#export_csv = result_df.to_csv (r'export_dataframe.csv', index = None, header=True) #Don't forget to add '.csv' at the end of the path
