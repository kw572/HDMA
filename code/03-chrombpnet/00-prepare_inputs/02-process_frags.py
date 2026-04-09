# Purpose: split fragments (grouped per cell type) into pseudoreplicates for peak calling
# Adapted from Salil Deshpande. We assume that fragments are already split by cell type.

import random
import sys

def load_reference_seqlevels(chromsizes_path):
    seqlevels = []
    with open(chromsizes_path, "r") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                seqlevels.append(parts[0])
    return set(seqlevels)


def normalize_chrom_name(chrom, reference_seqlevels):
    if chrom in reference_seqlevels:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom[3:]
    else:
        alt = f"chr{chrom}"
    if alt in reference_seqlevels:
        return alt
    return None

# parse the input file
filename = sys.argv[1]
cluster_sample = filename[:-4] if filename.endswith(".tsv") else filename
cluster = cluster_sample.split("__")[0] # get cluster
sample = cluster_sample.split("__")[-1] # get sample
chromsizes_path = sys.argv[3]
allowed_chros = load_reference_seqlevels(chromsizes_path)

print("--------------")
print(cluster_sample)
print(cluster)
print(sample)

# set up inputs and outputs
basedir = sys.argv[2]
frags_infile = f"{basedir}/fragments/{filename}"

out_pseudorep1_file = f"{basedir}/pseudorep1/{cluster}__{sample}.tsv"
out_pseudorep2_file = f"{basedir}/pseudorep2/{cluster}__{sample}.tsv"
out_pseudorepT_file = f"{basedir}/pseudorepT/{cluster}__{sample}.tsv"

# parse fragments files, splitting each one into pseudoreps
num_matches = 0
normalized_chroms = 0
with open(frags_infile, 'r') as f_f_in, open(out_pseudorep1_file, 'w') as f_p1_out, open(out_pseudorep2_file, 'w') as f_p2_out, open(out_pseudorepT_file, 'w') as f_pT_out:
	for line in f_f_in:
		chro, start, end, barcode = line.strip().split("\t")
		normalized = normalize_chrom_name(chro, allowed_chros)
		if normalized is None:
			continue
		if normalized != chro:
			normalized_chroms += 1
		chro = normalized
		start = int(start); end = int(end)
		if chro in allowed_chros:
			num_matches += 1
			# Output pseudorepT
			f_pT_out.write(f"{chro}\t{start}\t{start+1}\t{barcode}\n")
			f_pT_out.write(f"{chro}\t{end-1}\t{end}\t{barcode}\n")
			# Output start to p1/p2
			if random.random() < 0.5:
				f_p1_out.write(f"{chro}\t{start}\t{start+1}\t{barcode}\n")
			else:
				f_p2_out.write(f"{chro}\t{start}\t{start+1}\t{barcode}\n")
			# Output end to p1/p2
			if random.random() < 0.5:
				f_p1_out.write(f"{chro}\t{end-1}\t{end}\t{barcode}\n")
			else:
				f_p2_out.write(f"{chro}\t{end-1}\t{end}\t{barcode}\n")

print(f"{sample} {cluster}")
print(f"accepted_fragments={num_matches}")
print(f"normalized_chroms={normalized_chroms}")
