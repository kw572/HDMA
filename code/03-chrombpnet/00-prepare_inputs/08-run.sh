#!/bin/bash
set -euo pipefail

# Purpose: runner script which executes all commands to get GC-matched negative
# or background training regions, using
# an array job to limit the number of jobs in the queue.
# 1) 08-write_cmds.sh first loops over datasets and generates all commands
# 2) 08-run.sh is then submitted to the scheduler, which creates an array job over all cmds in 08-commands.sh
# 3) 08-get_negatives.sh is called by each command.

#SBATCH --job-name=08-get_negatives.sh
#SBATCH --output=../../logs/03-chrombpnet/00/08/%x-%j.out
#SBATCH --partition=main
#SBATCH --mem=20G
#SBATCH --array=1-1015:5
#SBATCH -n 5
#SBATCH --time=02:00:00

cmdfile="08-commands.sh"

if [[ ! -f "${cmdfile}" ]]; then
  echo "Command file not found: ${cmdfile}" >&2
  exit 1
fi

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  for i in {0..4}; do
    id=$((SLURM_ARRAY_TASK_ID+i))
    cmd=$(sed -n "${id}p" "${cmdfile}")
    [[ -n "${cmd}" ]] || continue
    bash -lc "${cmd}" &
  done
else
  while IFS= read -r cmd; do
    [[ -n "${cmd}" ]] || continue
    bash -lc "${cmd}"
  done < "${cmdfile}"
fi

# # important to make sure the job doesn't exit before the background tasks are done
# # https://www.sherlock.stanford.edu/docs/advanced-topics/job-management/#minimizing-the-number-of-jobs-in-queue
wait
