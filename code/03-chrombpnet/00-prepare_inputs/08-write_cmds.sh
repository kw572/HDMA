#!/usr/bin/bash

# Purpose: generate nonpeaks (GC-matched negatives), which will be used as background
# regions during model training. This script writes the commands, one per fold
# per cell type.
# Info: https://github.com/kundajelab/chrombpnet/wiki/Preprocessing#generate-non-peaks-background-regions

# SETUP ------------------------------------------------------------------------

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

export LC_ALL=C
export LC_CTYPE=C
export LANG=C

input_len=2114
stride=1000


# CONSTRUCT COMMANDS ------------------------------------------------

# get all clusters
datasets=$(
  find "${peaks_dir%/}" -maxdepth 1 -type f -name '*__peaks_overlap_filtered.narrowPeak' ! -name '._*' -print |
    while IFS= read -r path; do
      basename "${path}" __peaks_overlap_filtered.narrowPeak
    done
)

cmdfile=08-commands.sh
rm -f "${cmdfile}"
touch "${cmdfile}"

for cluster in ${datasets}; do
  
  echo "${cluster}"
    
  peak_file="${chrombpnet_peaks_dir}/${cluster}__peaks_bpnet.narrowPeak"
  if [[ ! -s "${peak_file}" ]]; then
    echo "${cluster} has no chrombpnet peaks; skipping"
    continue
  fi
  
  cluster_dir="${negatives_dir%/}/${cluster}"
  [[ -d "${cluster_dir}" ]] || mkdir -p "${cluster_dir}"
  
  for fold in {0..4}; do
    fold_name=fold_${fold}
    echo -e "\t${fold_name}"
  
    split_file="${split_dir%/}/${fold_name}.json"
  
    fold_dir="${cluster_dir%/}/${fold_name}"
    [[ -d "${fold_dir}" ]] || mkdir -p "${fold_dir}"
    fold_out="${fold_dir}/output"
  
    job_name="08-neg_${cluster}_${fold_name}"
    JOBSCRIPT="bash ./08-get_negatives.sh"
  
    # submit jobs
    done_file="${fold_dir}/output_negatives.bed"
    
    if [[ ! -f "$done_file" ]]; then
      rm -rf "${fold_dir}/output_auxiliary"
      printf "%s %q %q %q %q %q %q %q %q\n" \
        "${JOBSCRIPT}" "${ref_fasta}" "${chromsizes}" "${blacklist}" "${peak_file}" \
        "${input_len}" "${stride}" "${fold_out}" "${split_file}" >> "${cmdfile}"
    else
      echo "${cluster} ${fold} done; skipping"
    fi
    
  done

done
