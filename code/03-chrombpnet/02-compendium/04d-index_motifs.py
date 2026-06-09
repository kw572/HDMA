# Purpose: after identifying and annotating merged motifs
# we have certain patterns that we have manually identified as
# redundant. Here we assign an index to every unique motif.
# This operates on 04d-ChromBPNet_de_novo_motifs.tsv,
# and overwrites the same file.

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Assign a unique motif index per annotation label.")
parser.add_argument("--input", default="04d-ChromBPNet_de_novo_motifs.tsv")
args = parser.parse_args()

merged_anno = pd.read_csv(args.input, sep = "\t")
merged_anno['idx_uniq'] = merged_anno.groupby('annotation').ngroup() + 0
merged_anno.head()
print(len(set(merged_anno['idx_uniq'].to_list())))
merged_anno.to_csv(args.input, sep = "\t", index = False)
