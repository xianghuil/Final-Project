import pandas as pd
import numpy as np

if __name__ == '__main__':
	data = pd.read_csv('exp_seq.BRCA-KR.tsv', sep='\t')

	n_samples = data.shape[0]
	sample_ids = data['icgc_sample_id'].unique()
	gene_ids = data['gene_id'].unique()
	col_names = np.append('gene_id', sample_ids)
	cols = col_names.shape[0]
	rows = gene_ids.shape[0]

	# np.zeros((rows, cols))
	new_data = pd.DataFrame(np.zeros((rows, cols)), columns=col_names, index=gene_ids)
	new_data = new_data.astype('str')
	# print(new_data)

	for r in range(n_samples):
	    if r % 1000 == 0:
	        print(r)
	    col = data.iloc[r]['icgc_sample_id']
	    gene_id = data.iloc[r]['gene_id']
	#     print(gene_id)
	    new_data.loc[gene_id][col] = data.iloc[r]['normalized_read_count']
	    new_data.loc[gene_id]['gene_id'] = gene_id

	new_data.to_csv('data.csv', sep='\t', index=False)