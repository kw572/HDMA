#!/usr/bin/bash

set -euo pipefail

# source configuration variables
source ../config.sh

mkdir -p "${base_dir%/}/01-models/qc"

# count peaks
wc -l "${chrombpnet_peaks_dir}"/*.narrowPeak > "${base_dir%/}/01-models/qc/npeaks.tsv"
