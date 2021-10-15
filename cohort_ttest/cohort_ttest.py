import numpy as np
import pandas as pd
from scipy.stats import zscore, ttest_ind

# select genes
selected_genes = pd.read_excel('Selected_Pathway_Genes.xlsx',index_col=0,header=0)
# loading cohorts
voi_exprs = pd.read_excel('VoineaguRNAASD.xlsx',sheet_name='Expression_matrix')
voi_pheno = pd.read_excel('VoineaguRNAASD.xlsx',sheet_name='Phenotype')
voi_selected_genes = voi_exprs.loc[voi_exprs['Symbol'].isin(selected_genes['Genes']),:]
voi_selected_genes.dropna(inplace=True)
gse2mark_voi = dict(zip(voi_pheno['geo_accession'],voi_pheno['disease status:ch1']))
voi_ctrl_list = [key for key in gse2mark_voi if gse2mark_voi[key]=="controls"]
voi_asd_list = [key for key in gse2mark_voi if gse2mark_voi[key]=="autism"]
voi_selected_genes['Log2FC'] = np.log2(voi_selected_genes[voi_asd_list].mean(axis=1)/voi_selected_genes[voi_ctrl_list].mean(axis=1))
voi_selected_genes['t-stats'], voi_selected_genes['p-val'] = ttest_ind(voi_selected_genes[voi_asd_list],voi_selected_genes[voi_ctrl_list],axis=1)
voi_out = voi_selected_genes[["Symbol","Log2FC","p-val"]]

# Wright
wri_exprs = pd.read_excel('WrightRNAASD.xlsx'.format(datdir),sheet_name='Rpkm_matrix')
wri_pheno = pd.read_excel('WrightRNAASD.xlsx'.format(datdir),sheet_name='Phenotype')
wri_selected_genes = wri_exprs.loc[wri_exprs['hgnc_symbol'].isin(selected_genes['Genes']),:]
wri_selected_genes.dropna(inplace=True)
gse2mark_wri = dict(zip(wri_pheno['RiboRNum'].astype(str),wri_pheno['Dx']))
wri_ctrl_list = [str(key) for key in gse2mark_wri if gse2mark_wri[key]=="Control"]
wri_asd_list = [str(key) for key in gse2mark_wri if gse2mark_wri[key]=="Diseased"]
wri_selected_genes['Log2FC'] = np.log2(wri_selected_genes[wri_asd_list].mean(axis=1)/wri_selected_genes[wri_ctrl_list].mean(axis=1))
wri_selected_genes['t-stats'], wri_selected_genes['p-val'] = ttest_ind(wri_selected_genes[wri_asd_list],wri_selected_genes[wri_ctrl_list],axis=1)
wri_out = wri_selected_genes[["hgnc_symbol","Log2FC","p-val"]]

# Irimia
iri_exprs = pd.read_excel('IrimiaRNAASD.xlsx',sheet_name='Expression_matrix')
iri_pheno = pd.read_excel('IrimiaRNAASD.xlsx',sheet_name='Phenotype')
iri_sample = iri_pheno['title'].str.split('_',expand=True)
iri_sample['sample'] = iri_sample[1]+'_'+iri_sample[2].str.replace('-','.')+'_'+iri_sample[3]
iri_selected_genes = iri_exprs.loc[iri_exprs['hgnc_symbol'].isin(selected_genes['Genes']),:]
iri_selected_genes.dropna(inplace=True)
gse2mark_iri = dict(zip(iri_sample['sample'],iri_pheno['diagnosis:ch1']))
iri_ctrl_list = [key for key in gse2mark_iri if gse2mark_iri[key]=="CTL"]
iri_asd_list = [key for key in gse2mark_iri if gse2mark_iri[key]=="ASD"]
iri_selected_genes['Log2FC'] = np.log2(iri_selected_genes[iri_asd_list].mean(axis=1)/iri_selected_genes[iri_ctrl_list].mean(axis=1))
iri_selected_genes['t-stats'], iri_selected_genes['p-val'] = ttest_ind(iri_selected_genes[iri_asd_list],iri_selected_genes[iri_ctrl_list],axis=1)
iri_selected_genes['p-adj'] = p_adjust_bh(iri_selected_genes['p-val'])
iri_out = iri_selected_genes[["hgnc_symbol","Log2FC","p-val"]]

writer = pd.ExcelWriter("Cohorts_ttest.xlsx")
voi_out.to_excel(writer,sheet_name="Voi")
wri_out.to_excel(writer,sheet_name="Wri")
iri_out.to_excel(writer,sheet_name="Iri")
writer.save()
