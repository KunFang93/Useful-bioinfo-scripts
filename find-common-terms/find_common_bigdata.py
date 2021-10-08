import sys
import pandas as pd
from functools import reduce
from collections import Counter

print("Input file format: each row for a sample.txt file.")
try:
    names_list = sys.argv[1]
except IndexError:
    print("python find_commons.py genenames_list.txt")
    exit(1)

if len(sys.argv) >2:
    print("Input more than one files")
    exit(1)

# Reading
genenames_list = []
with open(names_list,'r') as files_list:
    for line in files_list:
        line_f = line.strip().split()
        if len(line_f) >1:
            print("Each row for one file!")
            exit(1)
        genenames_list.append(line_f[0])
files_list.close()

samples2id_dict = {}
for idx,sample in enumerate(genenames_list):
    samples2id_dict[sample] = idx

# Processing
genes_dict = {}
for idx,file in enumerate(genenames_list):
    print("Processing {}/{}".format(idx,len(genenames_list)),end="\r",flush=True)
    with open(file,'r') as tmp_f:
        for line in tmp_f:
            line_info = line.strip().split()
            line_gene = line_info[0]
            if len(line_info) >1:
                print("Each row for one gene!")
                exit(1)
            else:
                if line_gene not in genes_dict:
                    genes_dict[line_gene] = [0] * len(genenames_list)
                    genes_dict[line_gene][samples2id_dict[file]] = 1
                else:
                    genes_dict[line_gene][samples2id_dict[file]] = 1

df_merged = pd.DataFrame.from_dict(genes_dict,orient='index')
col_names = [sample.split('.')[0] for sample in genenames_list]
df_merged.columns = col_names
df_merged["Common Count"] = df_merged[col_names].sum(axis=1)

# Writing results
writer = pd.ExcelWriter("Common_genes.xlsx")
summary_dict = Counter(df_merged['Common Count'])
out_dict = {'Common Type':[],'Counts':[]}
for key in sorted(summary_dict.keys()):
    out_dict['Common Type'].append(">={}-Commons".format(key))
    current_df = df_merged[df_merged['Common Count'] >= key]
    out_dict['Counts'].append(len(current_df))
    current_df.to_excel(writer,'{}_commons'.format(key))
writer.save()
out_df = pd.DataFrame(out_dict)
out_df.to_excel("Common_genes_summary.xlsx",index=False)
