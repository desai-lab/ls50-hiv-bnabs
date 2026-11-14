import numpy as np 
import pandas as pd
import spoa
import sys
from pathlib import Path
from collections import Counter
import random
from rapidfuzz import process, fuzz
from rapidfuzz.distance import Levenshtein as RFLevenshtein
import networkx as nx
import argparse
import os
from numba import njit
import time

# Parameters
min_agree_vh = 0.75  # min agreement at every base to merge VH
min_agree_vl = 0.75  # min agreement at every base to merge VL
min_coverage = 10    # minimum coverage


@njit(cache=True)
def _sum_of_pairs_numba(msa_encoded, n_seqs, n_cols, match, mismatch, gap):
	total = 0
	gap_code = 4
	
	for col in range(n_cols):
		for i in range(n_seqs):
			c1 = msa_encoded[i, col]
			for j in range(i + 1, n_seqs):
				c2 = msa_encoded[j, col]
				if c1 == gap_code and c2 == gap_code:
					continue
				elif c1 == gap_code or c2 == gap_code:
					total += gap
				elif c1 == c2:
					total += match
				else:
					total += mismatch
	return total


def encode_msa_fast(msa):
	"""Fast MSA encoding using numpy vectorized operations."""
	arr = np.array([list(s) for s in msa], dtype='U1')
	encoded = np.full(arr.shape, 4, dtype=np.int8)  # Default to gap
	encoded[(arr == 'A') | (arr == 'a')] = 0
	encoded[(arr == 'C') | (arr == 'c')] = 1
	encoded[(arr == 'G') | (arr == 'g')] = 2
	encoded[(arr == 'T') | (arr == 't')] = 3
	return encoded


def sum_of_pairs_score(msa, match=1, mismatch=-1, gap=-2):
	if not msa:
		return 0
	encoded = encode_msa_fast(msa)
	return _sum_of_pairs_numba(encoded, len(msa), len(msa[0]), match, mismatch, gap)


def best_msa(sequences, n_trials=20, algorithm=1, **poa_kwargs):
	n_seqs = len(sequences)
	
	# Trivial cases - skip POA entirely
	if n_seqs == 0:
		return []
	if n_seqs == 1:
		return [sequences[0]]
	if len(set(sequences)) == 1:
		# All sequences identical
		return list(sequences)
	
	best_score = float('-inf')
	best_result = None
	best_order = None
	
	indices = list(range(n_seqs))
	
	# Adaptive trials: fewer for small groups
	effective_trials = min(n_trials, max(5, n_seqs * 2))
	
	recent_scores = []
	
	for trial in range(effective_trials):
		shuffled_idx = indices.copy()
		random.shuffle(shuffled_idx)
		shuffled_seqs = [sequences[i] for i in shuffled_idx]
		
		# Run POA alignment
		cons, msa = spoa.poa(shuffled_seqs, algorithm=algorithm, m=5, n=-4, g=-5, e=-3)
		
		# Fast scoring
		score = sum_of_pairs_score(msa)
		recent_scores.append(score)
		
		if score > best_score:
			best_score = score
			best_result = (cons, msa)
			best_order = shuffled_idx
		
		# Early termination if scores are identical
		if trial == 4:
			if len(set(recent_scores)) == 1:
				break
	
	cons, msa = best_result
	
	# Reorder MSA to match original input order
	inv_order = [0] * len(best_order)
	for new_pos, orig_idx in enumerate(best_order):
		inv_order[orig_idx] = new_pos
	
	msa_original_order = [msa[inv_order[i]] for i in range(n_seqs)]
	
	return msa_original_order


def check_aln(sequences_una,min_agree):
	sequences_una = list(sequences_una)
	if len(set(sequences_una)) == 1:
		return list(set(sequences_una))[0]
	
	sequences = best_msa(sequences_una)
	pos_msa = zip(*[list(sequences[i]) for i in range(len(sequences))])
	len_p = len(sequences)
	for p in pos_msa:
		if max(list(Counter(p).values())) / len_p < min_agree:
			return True
	return False

def fastq_qual_to_phred(qual_string):
	q = []
	for qual_string_i in qual_string:
		q.append(np.array([ord(c) - 33 for c in qual_string_i]))
	return q


def bayesian_consensus(seq_cc, qual_string, min_agree,
					   p_insertion=0.001, 
					   p_deletion=0.0005):
	"""
	Bayesian consensus using error probabilities.
	"""
	

	if len(set(seq_cc)) == 1:
		return list(set(seq_cc))[0]

	sequences = best_msa(seq_cc)

	seq_length = len(sequences[0])
	bases = ['A', 'C', 'G', 'T']
	base_to_idx = {b: i for i, b in enumerate(bases)}
	qualities = fastq_qual_to_phred(list(qual_string))
	
	consensus = []
	
	qual_indices = []
	for seq in sequences:
		idx_map = []
		original_pos = 0
		for char in seq:
			if char == '-':
				idx_map.append(None)
			else:
				idx_map.append(original_pos)
				original_pos += 1
		qual_indices.append(idx_map)
	
	for pos in range(seq_length):
		log_likes = np.zeros(5)
		
		for seq_idx, (seq, qual) in enumerate(zip(sequences, qualities)):
			obs_base = seq[pos].upper()
			is_gap = (obs_base == '-')
			
			if not is_gap:
				orig_pos = qual_indices[seq_idx][pos]
				if orig_pos is None or orig_pos >= len(qual):
					q = 20 # default, moderate quality
				else:
					q = qual[orig_pos]
				p_substitution = 10 ** (-q / 10)
				p_correct_base = 1 - p_substitution
			
			if is_gap:
				for true_idx in range(4):
					log_likes[true_idx] += np.log(p_deletion + 1e-10)
				log_likes[4] += np.log(1 - p_insertion + 1e-10)
			else:
				if obs_base not in base_to_idx:
					continue
				
				obs_idx = base_to_idx[obs_base]
				
				for true_idx in range(4):
					if true_idx == obs_idx:
						p_obs = (1 - p_deletion) * p_correct_base
					else:
						p_obs = (1 - p_deletion) * (p_substitution / 3)
					log_likes[true_idx] += np.log(p_obs + 1e-10)
				
				log_likes[4] += np.log(p_insertion * 0.25 + 1e-10)
		
		best_idx = np.argmax(log_likes)
		
		if best_idx < 4:
			consensus.append(bases[best_idx])
	
	return ''.join(consensus)


def check_n_clusters(sequences, qualities, max_mismatch=5):
	sequences = np.array(sequences)
	qualities = np.array(qualities)
	edges = [(i, i) for i in range(len(sequences))]
	for i1, s1 in enumerate(sequences):
		for i2, s2 in enumerate(sequences):
			if i1 > i2:
				if s1==s2:
					edges.append((i1, i2))
				elif RFLevenshtein.distance(s1, s2) <= max_mismatch:
					edges.append((i1, i2))

	G = nx.from_edgelist(edges)
	largest_cc = list(max(nx.connected_components(G), key=len))
	return len(largest_cc) / len(sequences), [sequences[i] for i in largest_cc], [qualities[i] for i in largest_cc], largest_cc


def process_chunk(df_chunk):
	"""Process a subset of super_bcs and return results."""
	superbcs = []
	purity = []
	vh_conc = []
	vl_conc = []
	vh_pur = []
	vl_pur = []
	bc_conc = []
	mini_conc = []
	coverage = []
	
	bcs_pass=[]
	bcs_pass_sbcs=[]
	c_notenoughcov = 0
	
	
	t_check=0
	t_lev=0
	t_bay=0

	for sbc, group in df_chunk.groupby(by='super_bc'):
		superbcs.append(sbc)
		coverage.append(len(group))
		
		if len(group) < min_coverage:
			purity.append(0)
			vh_conc.append(None)
			vl_conc.append(None)
			bc_conc.append(None)
			mini_conc.append(None)
			vh_pur.append(None)
			vl_pur.append(None)
			c_notenoughcov += 1
			continue

		start=time.time()
		cvh, cc_h, q_h, lcc = check_n_clusters(np.array(group['vh']), np.array(group['vh_qual']))
		
		if cvh < min_agree_vh:
			purity.append(0)
			vh_conc.append(None)
			vl_conc.append(None)
			bc_conc.append(None)
			mini_conc.append(None)
			vh_pur.append(None)
			vl_pur.append(None)
			continue

		cvl, cc_l, q_l, _ = check_n_clusters(np.array(group['vl']), np.array(group['vl_qual']))
		if cvl < min_agree_vl:
			purity.append(0)
			vh_conc.append(None)
			vl_conc.append(None)
			bc_conc.append(None)
			mini_conc.append(None)
			vh_pur.append(None)
			vl_pur.append(None)
			continue
		t_lev+=time.time()-start

		vh_pur.append(cvh)
		vl_pur.append(cvl)

		start=time.time()
		if check_aln(np.array(group['vh']),min_agree_vh):
			purity.append(0)
			vh_conc.append(None)
			vl_conc.append(None)
			bc_conc.append(None)
			mini_conc.append(None)
			continue
		
		if check_aln(np.array(group['vl']),min_agree_vl):
			purity.append(0)
			vh_conc.append(None)
			vl_conc.append(None)
			bc_conc.append(None)
			mini_conc.append(None)
			continue
		t_check+=time.time()-start

		

		start=time.time()
		vhc = bayesian_consensus(cc_h, q_h, min_agree_vh)

		
		if len(vhc) <= 1:
			purity.append(0)
			vh_conc.append(None)
			vl_conc.append(None)
			bc_conc.append(None)
			mini_conc.append(None)
			continue

		vlc = bayesian_consensus(cc_l, q_l, min_agree_vl)
		t_bay+=time.time()-start

		if len(vlc) <= 1:
			purity.append(0)
			vh_conc.append(None)
			vl_conc.append(None)
			bc_conc.append(None)
			mini_conc.append(None)
		else:
			vh_conc.append(vhc)
			vl_conc.append(vlc)

			# use the largest vh cc to compute the bc consensus
			bclii=np.array(group['bc'])
			bc_cc=[bclii[i] for i in lcc]
			bcliq=np.array(group['bc_qual'])
			bc_qual_cc=[bcliq[i] for i in lcc]
			bc_i = bayesian_consensus(bc_cc, bc_qual_cc, 0)

			bcs_pass+=bc_cc
			bcs_pass_sbcs+=[sbc]*len(bc_cc)

			bc_conc.append(bc_i)
			mini_conc.append(bc_i[0:7])
			purity.append(1)
	
	result_df = pd.DataFrame({
		'super_bc': superbcs,
		'vh_cons': vh_conc,
		'vl_cons': vl_conc,
		'bc_cons': bc_conc,
		'mini_bc': mini_conc,
		'purity': purity,
		'coverage': coverage,
		'vh_purity': vh_pur,
		'vl_purity': vl_pur
	})
	print(f'MSA check: {t_check}')
	print(f'Lev dist: {t_lev}')
	print(f'Bayesian: {t_bay}')

	bcs_pass_df=pd.DataFrame({
		'super_bc':bcs_pass_sbcs,
		'bc':bcs_pass
		})
	return result_df, c_notenoughcov, bcs_pass_df



def main():
	parser = argparse.ArgumentParser(description='Parallel barcode error correction')
	parser.add_argument('--task-id', type=int, required=True, 
						help='SLURM array task ID (1-indexed)')
	parser.add_argument('--num-tasks', type=int, required=True,
						help='Total number of parallel tasks')
	parser.add_argument('--output-dir', type=str, default='../chunks',
						help='Directory for chunk output files')
	args = parser.parse_args()
	
	# Create output directory if it doesn't exist
	os.makedirs(args.output_dir, exist_ok=True)
	
	print(f"Task {args.task_id} of {args.num_tasks} starting...")
	
	# Load data
	print("Loading data...")
	df2=pd.read_csv('../merged_clusters.csv')
	
	# Get all unique super_bcs
	all_super_bcs = np.array(pd.read_csv('../shuffled_sbcs.csv')['super_bc'])
	total_sbcs = len(all_super_bcs)
	print(f"Total super_bcs: {total_sbcs}")
	
	# Calculate chunk for this task (1-indexed task_id)
	chunk_size = (total_sbcs + args.num_tasks - 1) // args.num_tasks  # ceiling division
	start_idx = (args.task_id - 1) * chunk_size
	end_idx = min(start_idx + chunk_size, total_sbcs)
	
	my_super_bcs = all_super_bcs[start_idx:end_idx]
	print(f"Task {args.task_id}: Processing super_bcs {start_idx} to {end_idx} ({len(my_super_bcs)} super_bcs)")
	
	if len(my_super_bcs) == 0:
		print(f"Task {args.task_id}: No super_bcs to process, exiting.")
		# Write empty result file
		empty_df = pd.DataFrame(columns=['super_bc', 'vh_cons', 'vl_cons', 'bc_cons', 
										  'mini_bc', 'purity', 'coverage', 'vh_purity', 'vl_purity'])
		empty_df.to_csv(f"{args.output_dir}/chunk_{args.task_id}.csv", index=False)
		return
	
	# Process chunk
	df_chunk = df2[df2['super_bc'].isin(my_super_bcs)].dropna(subset=['bc','vh','vl'],how='any')
	del df2
	result_df, c_notenoughcov, bcs_pass_df = process_chunk(df_chunk)
	
	# Save chunk results
	output_file = f"{args.output_dir}/chunk_{args.task_id}.csv"
	result_df.to_csv(output_file, index=False)

	# Save chunk results--barcodes
	output_file = f"{args.output_dir}/bcs_{args.task_id}.csv"
	bcs_pass_df.drop_duplicates().to_csv(output_file, index=False)
	
	print(f"Task {args.task_id} complete:")
	print(f"  Processed: {len(my_super_bcs)} super_bcs")
	print(f"  Pure: {result_df['purity'].sum()}")
	print(f"  Not enough coverage: {c_notenoughcov}")
	print(f"  Output saved to: {output_file}")


if __name__ == '__main__':
	main()

