# Purpose: for a given cell type (cluster), load the GC-matched background
# regions and randomly select 100 regions per cell type. This creates a background
# set that can be reproducibly used for other analyses.
# Concatenate all the selected regions and save as a BED file.
#
# N.B.: The 100 backgrounds aren't evenly sampled across folds.
# Also, in the future, we should validate these backgrounds
# as being correctly one-hot-encodeable using tangermeme.utils._validate_input.
# Because when using these sequences downstream for in silico experiments,
# sequences will be validated that way. Jacob's suggestion: load all loci
# and then subset to valid ones without missing characters:
# X[X.sum(dim=(1, 2)) == X.shape[-1]]

import os
import sys
import pandas as pd
import numpy as np
import argparse

def parse_args():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--cluster", type=str, default=None, help="Cluster for which to select background regions.")
	parser.add_argument("--negatives-dir", type=str, default=None, help="Base directory for the negatives, expected to contain a folder matching --cluster.")
	parser.add_argument("--output-dir", type=str, default=None, help="Output directory for the selected background regions.")
	parser.add_argument("--num-folds", type=int, default=5, help="Number of folds to sample from.")
	
	args = parser.parse_args()

	print(args)

	return args


def main(args):

	cluster = args.cluster
	negatives_dir = args.negatives_dir
	output_dir = args.output_dir
	folds = list(range(args.num_folds))
	
	print("------------" + cluster + "------------")

	nonpeaks = []

	# load GC-matched background sequences
	for fold in folds:
		
		nonpeaks_bed = os.path.join(negatives_dir, f"{cluster}/fold_{fold}/output_negatives.bed")
		
		# read in BED file
		nonpeaks_fold = pd.read_csv(nonpeaks_bed, sep="\t", header=None)
		
		# select 20 random indices
		n_to_sample = min(20, len(nonpeaks_fold.index))
		if n_to_sample == 0:
			continue
		indices = np.random.choice(nonpeaks_fold.index, n_to_sample, replace=False)

		# subset the peaks to the indices
		nonpeaks_subset = nonpeaks_fold.iloc[indices]

		# append
		nonpeaks.append(nonpeaks_subset)

	# concatenate
	if len(nonpeaks) == 0:
		raise ValueError(f"No negatives available to subset for cluster {cluster}")

	nonpeaks = pd.concat(nonpeaks)

	# write to file
	nonpeaks.to_csv(os.path.join(output_dir, f"{cluster}__output_negatives_rand100.bed"), sep="\t", header=False, index=False)


if __name__ == "__main__":
	args = parse_args()
	main(args)
