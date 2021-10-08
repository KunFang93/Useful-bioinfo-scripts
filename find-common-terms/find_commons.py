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

# Processing
genes_binary_df_list = []
colnames = []
for file in genenames_list:
    print("reading {}".format(file))
    col_name = file.split('/')[-1].split('.')[0]
    colnames.append(col_name)
    current_df = pd.read_table(file,sep='\t',header=None)
    current_df[col_name] = 1
    current_df.columns = ['Genes',col_name]
    genes_binary_df_list.append(current_df)

df_merged = reduce(lambda left,right: pd.merge(left,right,on=['Genes'],how='outer'), genes_binary_df_list).fillna(0)
df_merged["Common Count"] = df_merged[colnames].sum(axis=1)

# Writing results
writer = pd.ExcelWriter("Common_genes.xlsx")
summary_dict = Counter(df_merged['Common Count'])
out_dict = {'Common Type':[],'Counts':[]}
for key in sorted(summary_dict.keys()):
    out_dict['Common Type'].append("{}-Commons".format(key))
    current_df = df_merged[df_merged['Common Count'] >= key]
    out_dict['Counts'].append(len(current_df))
    current_df.to_excel(writer,'{}_commons'.format(key))
writer.save()
out_df = pd.DataFrame(out_dict)
out_df.to_excel("Common_genes_summary.xlsx",index=False)
