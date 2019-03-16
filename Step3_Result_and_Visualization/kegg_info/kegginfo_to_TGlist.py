import pandas as pd
import numpy as np
import re
import argparse

parser = argparse.ArgumentParser(description='Adding Kegg_information to the given TGlist, kegg_ath_gene.txt must be in the same directory')
parser.add_argument('TGlist', type=str, help='0.AT1G72740.TG.cold.D2.1')
parser.add_argument('-o','--output', type=str, help='output_file_name')
args = parser.parse_args()

# importing kegg_info and TG_list
kegg_info = pd.read_csv('../kegg_ath_gene.txt',header=None, delimiter='\t')
TG_list = pd.read_csv(str(args.TGlist), delimiter='\t') 

# striping unnecessary string in the kegg_info's gene_name column
for g in range(kegg_info.shape[0]):
	kegg_info.iloc[g,0] = kegg_info.iloc[g,0].strip('ath:')

# turning kegg_info into a dictionary
index = kegg_info.iloc[:,0].isin(TG_list.iloc[:,0])
filtered_kegg_info = kegg_info[index]
				    
filtered_kegg_info.columns = ['key','value']

kegg_dict = {k: g['value'].apply(lambda x: x.split(';')[0]) for k,g in filtered_kegg_info.groupby('key')}
print(kegg_dict)

# adding gene information to the corresponding TGs in the TG_list as the 2rd column
info_col = pd.Series(np.zeros(TG_list.shape[0]))

for i in range(TG_list.shape[0]):
	a = str(kegg_dict[str(TG_list.iloc[i,0])]) # kegg_information_part -> a
	info_col[i] = re.sub('[0-9]', '', a).strip().rstrip('\nName: value, dtype: object')
	
TG_list['info'] = info_col
	
TG_list.to_csv(''.join([str(args.TGlist),'_w_description']),sep='\t',header=True, index=False)

