def plot_heterogeneity():
	import pandas as pd
	import numpy as np
	import glob, os
	import csv
	import math
	import matplotlib
	import matplotlib.pyplot as plt
	import seaborn as sns
	import copy
	from tqdm import tqdm_notebook
	from itertools import tee
	import re
	path = "outputs"
	csv_files = []
	for file in os.listdir(path):
	    if file.endswith(".csv"):
	        csv_files.append(os.path.join(path, file))
	csv_files = sorted(csv_files)

	### SET MOLECULAR WEIGHTS OF MONOMERS IN ORDER OF LABEL (labels are 1 indexed)
	Mw = [100.121,500,198.3,246.32, 104.15]
	HLBs = [8.45, 11.42, 5.125, 18.5, 4.865]  # 5th one is made up for Styrene - NEED TO REPLACE.
	#Read csv output files of 'n' batches and convert to n-element pd dataframe
	dfs = []
	seq_lens = []

	for m in range(len(csv_files)):
	    seq = []
	    row_size = [] 
	    
	    with open(csv_files[m], mode = 'r') as file:
	        csvFile = csv.reader(file)
	        i = 0
	        for lines in csvFile:
	            i = i + 1
	            seq.append(lines)
	            row_size.append(len(lines))
	        df = pd.DataFrame(seq)
	        dfs.append(df)
	        seq_lens.append(row_size)
	### CALCULATE AVERAGE DP, Mn, Mw & PDI OF ALL SEQUENCES

	names = []
	num_seqs = []
	avg_Mns = []
	avg_Mws = []
	PDIs = []
	avg_DPs = []
	seqs = []

	for m in tqdm_notebook(range(len(csv_files))):
	    seq = []

	    with open(csv_files[m], mode = 'r') as file:
	        csvFile = csv.reader(file)
	        i = 0
	        for lines in csvFile:
	            i = i + 1
	            seq.append(lines)
	    raw_seq = copy.deepcopy(seq)
	    seqs.append(raw_seq)
	    seqlen = np.zeros(i)
	    seqweight = np.zeros(i)
	    seqweight2 = np.zeros(i)

	    for j in range(0,i):
	        seqlen[j] = len(seq[j])
	        for k in range(len(seq[j])):
	            seq[j][k] = Mw[int(seq[j][k])-1]
	        seqweight[j] = sum(seq[j])
	        seqweight2[j] = (sum(seq[j])**2)

	    AvgMw = sum(seqweight2)/sum(seqweight)
	    AvgMn = sum(seqweight)/len(seqweight)
	    PDI = AvgMw/AvgMn
	    DP = sum(seqlen)/len(seqlen)
	    
	    names.append(csv_files[m].replace(path, '').replace('/', '').replace('.csv', ''))
	    n_mons = []
	    MR1s = []
	    MR2s = []
	    MR3s = []
	    MR4s = []
	    MR5s = []
	    for name in names:
	        try:
	            # n_mons, MR1, ... MRN, n_chains, DP, conv, CTP, Filt
	            num_vals = re.findall(r'\d+', name)
	            n_mons.append(int(num_vals[0]))
	            MR1s.append(int(num_vals[1]))
	            MR2s.append(int(num_vals[2]))
	            if len(num_vals) > 8:
	                MR3s.append(int(num_vals[3]))
	            else:
	                MR3s.append(0)
	            if len(num_vals) > 9:
	                MR4s.append(int(num_vals[4]))
	            else:
	                MR4s.append(0)
	            if len(num_vals) > 10:
	                MR5s.append(int(num_vals[5]))
	            else:
	                MR5s.append(0)
	        except IndexError as err:
	            n_mons.append(4)
	            MR1s.append(0)
	            MR2s.append(0)
	            MR3s.append(0)
	            MR4s.append(0)
	            MR5s.append(0)
	    num_seqs.append(i)
	    avg_Mns.append(AvgMn)
	    avg_Mws.append(AvgMw)
	    PDIs.append(PDI)
	    avg_DPs.append(DP)
	    
	d = {'Batch Name': names, 'Mons': n_mons, 'MR1': MR1s, 'MR2': MR2s, 'MR3': MR3s, 'MR4': MR4s, 'MR5': MR5s, 'NumSeqs': num_seqs, 'Avg Mn': avg_Mns, 'Avg Mw': avg_Mws, 'PDI': PDIs, 'Avg DP': avg_DPs}
	df_prop = pd.DataFrame(data=d)
	### VISUALIZE SEQUENCE SIMULATED
	batch_no = 0
	seq_no = 0
	num_seqs = 10

	chains = []
	for i in range(num_seqs):
	    chains.append([int(x) if x != None else 0 for x in dfs[batch_no].iloc[seq_no + i]])

	plt.figure(num=None, figsize=(60, 2), dpi=80, facecolor='w', edgecolor='k')
	ax = sns.heatmap(chains, vmin=0, vmax=df_prop.iloc[batch_no]['Mons'], linewidth=0.2, xticklabels=False, yticklabels=False)
	plt.show()