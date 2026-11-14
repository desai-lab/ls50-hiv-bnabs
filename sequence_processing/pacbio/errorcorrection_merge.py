import pandas as pd
import numpy as np
import argparse
import glob
import os

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

def find_chunks(args,filename):
	# Find all chunk files
	chunk_files = sorted(glob.glob(f"{args.chunk_dir}/{filename}_*.csv"))
	print(f"Found {len(chunk_files)} chunk files")
	
	# Verify we have all expected chunks
	expected_chunks = set(range(1, args.num_tasks + 1))
	found_chunks = set()
	for f in chunk_files:
		chunk_num = int(os.path.basename(f).replace(filename+'_', '').replace('.csv', ''))
		found_chunks.add(chunk_num)
	
	missing_chunks = expected_chunks - found_chunks
	if missing_chunks:
		print(f"WARNING: Missing chunks: {sorted(missing_chunks)}")
		print("Some tasks may have failed. Proceeding with available chunks...")
	
	# Load and concatenate all chunks
	dfs = []
	for f in chunk_files:
		df = pd.read_csv(f)
		dfs.append(df)
		print(f"  Loaded {f}: {len(df)} rows")
	
	if not dfs:
		print("ERROR: No chunk files found!")
		return
	
	return pd.concat(dfs, ignore_index=True)

def main():
	parser = argparse.ArgumentParser(description='Merge parallel error correction results')
	parser.add_argument('--chunk-dir', type=str, default='../chunks',
						help='Directory containing chunk output files')
	parser.add_argument('--num-tasks', type=int, required=True,
						help='Total number of parallel tasks (to verify all completed)')
	args = parser.parse_args()
	
	print("Merging chunk results...")
	
	pur_df=find_chunks(args,'chunk')
	
	print(f"\nCombined results: {len(pur_df)} total super_bcs")
	
	# Calculate statistics
	total_pure = int(pur_df['purity'].sum())
	total_impure = int((1 - pur_df['purity']).sum())
	
	print(f"\nTotal number of super_bcs remaining: {total_pure}")
	print(f"Total number of super_bcs thrown out: {total_impure}")
	
	# Load original ec_bcs for final output
	ec_bcs = parse_starcode_output('../clustered_barcodes.txt')
	
	# Merge and save outputs
	ec_bcs2 = ec_bcs.merge(pur_df[['super_bc', 'purity']], on='super_bc', how='left')
	ec_bcs2.to_csv('../superbarcodes_all.csv', index=False)
	print(f"Saved: ../superbarcodes_all.csv ({len(ec_bcs2)} rows)")
	
	ec_bcs2[ec_bcs2['purity'] == 1].to_csv('../superbarcodes_pure.csv', index=False)
	print(f"Saved: ../superbarcodes_pure.csv ({len(ec_bcs2[ec_bcs2['purity'] == 1])} rows)")
	
	pur_df.dropna(subset='bc_cons').to_csv('../consensus_seqs.csv', index=False)
	print(f"Saved: ../consensus_seqs.csv ({len(pur_df.dropna(subset='bc_cons'))} rows)")
	
	bc_df=find_chunks(args,'bcs').drop_duplicates()

	bc_df.to_csv('../barcodes_pure.csv',index=False)

	print('\nFinished!')


if __name__ == '__main__':
	main()

