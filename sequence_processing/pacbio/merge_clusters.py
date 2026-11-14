import numpy as np 
import pandas as pd
import random

def parse_starcode_output(starcode_file):
	clusters = []
	
	with open(starcode_file, 'r') as f:
		for i, line in enumerate(f):
			parts = line.strip().split('\t')
			canonical = parts[0]
			count = int(parts[1])
			
			if len(parts) > 2:
				members = parts[2].split(',')
			else:
				members = [canonical]
			
			super_bc = f'BC_{i+1}'
			
			for member in members:
				clusters.append({
					'super_bc': super_bc,
					'bc': member,
					'cluster_size': count
				})
	
	return pd.DataFrame(clusters)


ec_bcs = parse_starcode_output('../clustered_barcodes.txt')
df = pd.read_csv('../grab_barcodes_filtered.csv').dropna(subset='bc').astype(str)
df2 = df.merge(ec_bcs, on='bc', how='left')

df2.to_csv('../merged_clusters.csv',index=False)

all_super_bcs = list(df2['super_bc'].dropna().unique())
random.shuffle(all_super_bcs)

with open('../shuffled_sbcs.csv', 'w') as file:
	file.write(f"super_bc\n")
	for item in all_super_bcs:
		# Write the item followed by a newline character
		file.write(f"{item}\n")