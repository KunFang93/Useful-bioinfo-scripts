import sys
import pybedtools as pybed
import pandas as pd

try:
    meta_info_f = sys.argv[1]
    expand_len = int(sys.argv[2])
    cutoff = float(sys.argv[3])
except IndexError:
    print("Usage: python find_chip_exon.py meta_info.xlsx exon_expand_length(INT) cutoff(FLOAT,[0,1])")
    exit(1)

chip_suffix = '_chip.bed'
exon_suffix = '_exon.bed'
datdir= './'
info_dict = {}
meta_info = pd.read_excel('{}/{}'.format(datdir,meta_info_f))
bimat_col = meta_info['Chip'] + '_' + meta_info['RNA']
total_len = len(meta_info)
for idx,value in meta_info.iterrows():
    print("{}/{} comparisons".format(idx+1,total_len),flush=True,end='\r')
    current_chip = '{}{}'.format(value['Chip'],chip_suffix)
    current_exon = '{}{}'.format(value['RNA'],exon_suffix)
    current_exon_df = pd.read_csv('{}/{}'.format(datdir,current_exon),sep='\t',header=0)
    current_exon_df[['start']],current_exon_df[['end']] = current_exon_df[['start']] - expand_len, current_exon_df[['end']] + expand_len
    current_chip_df = pd.read_csv('{}/{}'.format(datdir,current_chip),sep='\t',header=0)
    current_exon_bedtool = pybed.BedTool.from_dataframe(current_exon_df)
    current_chip_bedtool = pybed.BedTool.from_dataframe(current_chip_df)
    current_inter = current_exon_bedtool.intersect(current_chip_bedtool,wa=True,f=cutoff)
    current_inter_df = current_inter.to_dataframe()
    current_inter_df.columns = current_exon_df.columns
    for idy,value_y in current_inter_df.iterrows():
        current_exon_name = value_y['exon']
        if current_exon_name in info_dict:
            info_dict[current_exon_name][idx] = 1
        else:
            info_dict[current_exon_name] = [0]*total_len
            info_dict[current_exon_name][idx] = 1

info_df = pd.DataFrame.from_dict(info_dict,orient='index',columns=bimat_col)
info_df.to_excel('chip_exon_overlap_info.xlsx')
