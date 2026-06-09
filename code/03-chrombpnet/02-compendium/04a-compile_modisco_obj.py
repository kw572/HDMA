# Purpose: Compile merged motifs
# This script combines several h5 modisco objects into one.
# Specifically, we use the h5 modisco files created from merging modisco
# motifs within gimme cluster patterns of motifs.
# https://github.com/austintwang/finemo_gpu/blob/main/src/finemo/data_io.py#L168
# NOTE: this is relatively fast. It takes about 1-2 minutes to compile ~800 patterns
# from ~100 objects.


import argparse
import h5py as h5
import numpy as np
import pandas as pd
import os
import modiscolite.report


def parse_args():
    parser = argparse.ArgumentParser(description="Compile merged MoDISco objects into one H5/TSV bundle.")
    parser.add_argument("--modisco-merged-dir", required=True, type=str)
    parser.add_argument("--cluster-key", required=True, type=str)
    parser.add_argument("--out-h5", required=True, type=str)
    parser.add_argument("--out-tsv", required=True, type=str)
    return parser.parse_args()


args = parse_args()

compiled_modisco_h5_path = args.out_h5
compiled_modisco_tsv_path = args.out_tsv
cluster_key_path = args.cluster_key

cluster_key = pd.read_csv(cluster_key_path, sep="\t", header=None, names=["gimme_cluster", "pattern_class", "component_patterns", "batch", "n_component_patterns", "n_seqlets"])


# COMPILE PATTERNS ----------------------------------------------
# create a dict from pattern to modisco merged file
input_files = {os.path.basename(root): os.path.join(root, file) for root, dirs, files in os.walk(args.modisco_merged_dir)
               for file in files if file == "merged_modisco.h5"}

# check we have the same number of input files as clusters in the cluster key
len(input_files)

# open the compiled modisco h5
with h5.File(compiled_modisco_h5_path, "w") as compiled_modisco:

    # set up groups in h5 file
    h5_pattern_groups = {
        "pos_patterns": compiled_modisco.create_group('pos_patterns'),
        "neg_patterns": compiled_modisco.create_group('neg_patterns')
    }
    
    # loop over the input files containing merged patterns derived from collapsing/merging similar patterns
    # within each gimme cluster
    for pattern_cluster in input_files:
        
        print("@ processing:", pattern_cluster, " | input file:", input_files[pattern_cluster])

        # an input modisco obj should only have one type of pattern class, either pos or neg patterns.
        with h5.File(input_files[pattern_cluster]) as modisco_obj:

            # figure out what class of patterns the object contains
            if "pos_patterns" in modisco_obj:
                pattern_class = "pos_patterns"
            elif "neg_patterns" in modisco_obj:
                pattern_class = "neg_patterns"

            current_group = h5_pattern_groups[pattern_class]

            # iterate over patterns in the obj and add them to the compiled object
            for pattern in modisco_obj[pattern_class].keys():
                
                current_group.create_dataset(f"{pattern}/contrib_scores",
                                             data=modisco_obj[pattern_class][pattern]["contrib_scores"])
                current_group.create_dataset(f"{pattern}/sequence",
                                             data=modisco_obj[pattern_class][pattern]["sequence"])
                current_group.create_dataset(f"{pattern}/hypothetical_contribs",
                                             data=modisco_obj[pattern_class][pattern]["hypothetical_contribs"])
                # add in placeholder value for number of seqlets
                current_group.create_dataset(f"{pattern}/seqlets/n_seqlets", data=np.array([1]))


# also combine all the TSVs into one
merged_report_long = pd.concat([pd.read_csv(os.path.join(root, file), sep = "\t")
                                for root, dirs, files in os.walk(args.modisco_merged_dir)
               for file in files if file == "merged_modisco.tsv"])

merged_report_long.to_csv(compiled_modisco_tsv_path, sep = "\t", index = False)

print("@ done.")
