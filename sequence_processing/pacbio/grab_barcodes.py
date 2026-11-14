import numpy as np 
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import edlib

def fuzzy_find(pattern, text, max_errors=1):
    result = edlib.align(
        pattern, 
        text, 
        mode="HW",  # "HW" = infix/substring matching
        task="locations",  # we want positions, not just score
        k=max_errors  # max edit distance allowed
    )
    
    
    return result["editDistance"], result["locations"]

def pull_seq_bw(seq,qual,upstream,downstream,l):
	loc1=seq.find(upstream)+l
	if loc1==-1+l:
		ed,loc1v=fuzzy_find(upstream, seq, max_errors=1)
		if ed>1 or ed<0:
			return None, None
		loc1=loc1v[0][1]+1

	loc2=seq.find(downstream)
	if loc2==-1:
		ed,loc2v=fuzzy_find(downstream, seq, max_errors=1)
		if ed>1 or ed<0:
			return None, None
		loc2=loc2v[0][0]
	return seq[loc1:loc2], qual[loc1:loc2]

fastq_gen = FastqGeneralIterator(open('../pacbio-filtered.fastq'))

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
		





