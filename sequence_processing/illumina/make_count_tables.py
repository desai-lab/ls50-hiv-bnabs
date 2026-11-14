import numpy as np 
import pandas as pd


########## CHANGE ###############
metadata=pd.read_csv('../seq_metadata.csv')
input_dir='/n/netscratch/desai_lab/Lab/joao/ls50/umi_counts/'
input_dir2='/n/netscratch/desai_lab/Lab/joao/ls50/umi_counts_subsample/'
output_dir='/n/holylabs/desai_lab/Lab/joao/ls50/illumina/'

#experiments=list(metadata['experiment'].unique())



for experiment,group1 in metadata.groupby(by='experiment'):
	print(group1)
	c=0
	for j,row in group1.iterrows():
		try:
			df0=pd.read_csv(input_dir+str(row['name'])+'_umi_counts2.csv').rename(columns={'umi_counts':str(row['name'])})
		except:
			df0=pd.read_csv(input_dir2+str(row['name'])+'_umi_counts2.csv').rename(columns={'umi_counts':str(row['name'])})
		if c==0:
			df=df0
		else:
			df=df.merge(df0,on='super_bc',how='outer')
		c+=1
	
	df.fillna(0).to_csv(output_dir+experiment+'_counts.csv',index=False)
	del df
