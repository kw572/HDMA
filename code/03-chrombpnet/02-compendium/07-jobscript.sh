#!/usr/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/07/%x-%j.out
#SBATCH --partition=main
#SBATCH -t 96:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=6G

# get args
dataset=${1}
out=${2}
fragments=${3}
peaks=${4}
fasta=${5}
hits=${6}

set -euo pipefail

# setup environment
NUCLEOATAC_VENV="${NUCLEOATAC_VENV:-/scratch1/${USER}/venvs/nucleoatac310}"
NUCLEOATAC_PREFIX="${NUCLEOATAC_PREFIX:-/scratch1/${USER}/envs/nucleoatac27}"
if [[ ! -x "${NUCLEOATAC_VENV}/bin/activate" && -x "/scratch1/${USER}/venvs/nucleoatac/bin/activate" ]]; then
  NUCLEOATAC_VENV="/scratch1/${USER}/venvs/nucleoatac"
elif [[ ! -x "${NUCLEOATAC_VENV}/bin/activate" && -x "$HOME/venvs/nucleoatac/bin/activate" ]]; then
  NUCLEOATAC_VENV="$HOME/venvs/nucleoatac"
fi
if [[ -x "${NUCLEOATAC_VENV}/bin/activate" ]]; then
  # shellcheck disable=SC1090
  source "${NUCLEOATAC_VENV}/bin/activate"
elif [[ -x "${NUCLEOATAC_PREFIX}/bin/nucleoatac" ]]; then
  export PATH="${NUCLEOATAC_PREFIX}/bin:${PATH}"
else
  module load python/2.7.13
  eval "$(conda shell.bash hook)"
  conda activate nucleoatac
fi

# check version
nucleoatac --version

echo $out

cores=12

# run NucleoATAC


# ------------------------------------------------------------------------------
if [[ -f "${out}.occpeaks.bed.gz" ]]; then

  echo -e "@ Found OCC output"
  
else
  
  nucleoatac occ \
  	--fragments $fragments \
	--bed $peaks \
	--fasta $fasta \
	--out $out \
	--cores $cores \
	--upper 300

  nucleoatac vprocess --sizes ${out}.nuc_dist.txt \
  	--out $out \
  	--upper 300
	
fi


# ------------------------------------------------------------------------------
if [[ -f "${out}.nucpos.bed.gz" ]]; then

  echo -e "@ Found NUC output"
  
else
  
  nucleoatac nuc \
  	--fragments $fragments \
  	--bed $peaks \
  	--vmat ${out}.VMat \
  	--out $out \
  	--fasta $fasta \
  	--sizes ${out}.fragmentsizes.txt \
  	--occ_track ${out}.occ.bedgraph.gz \
  	--fasta $fasta \
  	--cores $cores
  	
fi


# ------------------------------------------------------------------------------
if [[ -f "${out}.nucmap_combined.bed.gz" ]]; then

  echo -e "@ Found MERGE output"
  
else
  
  nucleoatac merge \
  	--out $out \
  	--occpeaks ${out}.occpeaks.bed.gz \
  	--nucpos ${out}.nucpos.bed.gz
  	
fi
	

# ------------------------------------------------------------------------------
python 07-calc_dyad_dist.py ${out}.occpeaks.bed.gz ${out} --max_dist 250 --k_nearest 100


# ------------------------------------------------------------------------------
# annotate hits by distance to dyads and plot as heatmap
export R_LIBS_USER=$RENV_SJ
Rscript 07-plot_nucleoatac.R "${out}.occpeaks.bed.gz" ${hits} ${out} "${out}.dyad_binned_distances.tsv"
	
echo "@ done."
