import numpy as np 
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import edlib
from rapidfuzz import process, fuzz
from rapidfuzz.distance import Levenshtein
import gzip
import sys

max_ind_dist=1 # maximum distance of observed inline index to expected sequence


def fuzzy_find(pattern, text, max_errors=1):
	result = edlib.align(
		pattern, 
		text, 
		mode="HW",  # "HW" = infix/substring matching
		task="locations",  # we want positions, not just score
		k=max_errors  # max edit distance allowed
	)
	
	return result["editDistance"], result["locations"]

def pull_seq_bw(seq,upstream,downstream,l):
	loc1=seq.find(upstream)+l
	if loc1==-1+l:
		ed,loc1v=fuzzy_find(upstream, seq, max_errors=1)
		if ed>1 or ed<0:
			return None
		loc1=loc1v[0][1]+1

	loc2=seq.find(downstream)
	if loc2==-1:
		ed,loc2v=fuzzy_find(downstream, seq, max_errors=1)
		if ed>1 or ed<0:
			return None
		loc2=loc2v[0][0]
	return seq[loc1:loc2]

# INPUT
filename=sys.argv[1]#'0_JAsce_LA_6_3_S7_L003'
filedir=sys.argv[2]
fileloc=filedir+filename
outdir=sys.argv[3]
# CHANGE DIR
metadata=pd.read_csv('../seq_metadata.csv').set_index('name')


name0=fileloc.split('/')[-1]

name='_'.join(name0.split('0_JAsce_')[1].split('_')[:-2])

fastq_gen = FastqGeneralIterator(open(fileloc))
contam_seq='AGGCCCG'

upstream='TGGCGTAATAGCGAAG'
downstream='CTGAATGGCGAATGGACGC'
upstream_l=len(upstream)

reads_not_matching_index=0
reads_contam=0

bc_counts={}
rest_of_reads=0
total_reads=0




#print(metadata.loc[name])

r1_ind=str(metadata.loc[name]['r1 index'])
r2_ind=str(metadata.loc[name]['r2 index'])

r1_ind_l=len(r1_ind)
r2_ind_l=len(r2_ind)


for entry in fastq_gen:
	total_reads+=1
	# check inline indices
	seq=str(entry[1])
	r1_ind_obs=seq[8:r1_ind_l+8]
	
	if r1_ind_obs==r1_ind:
		pass
	elif Levenshtein.distance(r1_ind_obs,r1_ind)>max_ind_dist:
		reads_not_matching_index+=1
		continue

	r2_ind_obs=seq[-8-r2_ind_l:-8]
	if r2_ind_obs==r2_ind:
		pass
	elif Levenshtein.distance(r2_ind_obs,r2_ind)>max_ind_dist:
		reads_not_matching_index+=1
		continue
	

	umi1=seq[:8]
	umi2=seq[-8:]
	bc=pull_seq_bw(seq,upstream,downstream,upstream_l)
	if bc==None:
		reads_not_matching_index+=1
		continue

	minibc=bc[:7]
	if Levenshtein.distance(minibc,contam_seq)<=1:
		#print(bc)
		reads_contam+=1
		continue

	new_key=(umi1,umi2,bc)
	if new_key in bc_counts:
		bc_counts[new_key]+=1
	else:
		bc_counts[new_key]=1
	rest_of_reads+=1

print('Total Reads:',total_reads)
print('Passed Reads:',rest_of_reads)
print('Fraction Pass:',rest_of_reads/total_reads)
print('Contaminated Reads:',reads_contam)
print('Reads with wrong index:',reads_not_matching_index)


with open(f"{outdir}{name}_init_counts.csv", "w") as f:
	f.write('umi1,umi2,bc,counts\n')
	for key in bc_counts:
		f.write(f'{key[0]},{key[1]},{key[2]},{bc_counts[key]}\n')

'''
upstream_bc='TGGCGTAATAGCGAAG'
downstream_bc='CTGAATGGCGAAT'
upstream_VH='GTGGTGGTGGTTCTGCTAGC'
downstream_VH='AGTGGTGGAGGCGGATCAGG'
upstream_VL='GAAGTGGCGGCGGAGGTTCA'
downstream_VL='GAACAAAAGCTTATTTCTGA'

lbc=len(upstream_bc)
lvh=len(upstream_VH)
lvl=len(upstream_VL)
with open("../grab_barcodes.csv", "w") as f:
	f.write('id,bc,bc_qual,vh,vh_qual,vl,vl_qual\n')
	for entry in fastq_gen:
		bc,bc_qual = pull_seq_bw(entry[1],entry[2],upstream_bc,downstream_bc,lbc)
		if bc==None:
			continue
		vh,vh_qual = pull_seq_bw(entry[1],entry[2],upstream_VH,downstream_VH,lvh)
		if vh==None:
			continue
		vl,vl_qual = pull_seq_bw(entry[1],entry[2],upstream_VL,downstream_VL,lvl)
		if vl==None:
			continue
		#write all strings + identity to a file
		f.write(entry[0]+','+bc+','+bc_qual+','+vh+','+vh_qual+','+vl+','+vl_qual+'\n')
		
'''



