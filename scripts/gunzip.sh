#!/usr/bin/env bash
#SBATCH -J gunzip         ## $SLURM_JOB_NAME
#SBATCH -p compute        ## which set of physical compute servers to use
#SBATCH -q batch          ## quality of service (set of resource limits)
#SBATCH -t 72:00:00       ## time limit (max is 72h for '-q batch')
#SBATCH -c 1              ## numbers of CPUs (actually threads/vCPUs)
#SBATCH --mem 1g          ## RAM; 1 gigabyte
#SBATCH -o %x.%j.out      ## redirect stdout to $SLURM_JOB_NAME.$SLURM_JOB_ID.out
#SBATCH -e %x.%j.err      ## redirect stderr to $SLURM_JOB_NAME.$SLURM_JOB_ID.err

## $# (number of arguments) equals 1 or print usage to stderr (>&2) and exit non-0:

[ $# -eq 1 ] || {  
  echo "Usage: sbatch gunzip.sh <input.gz>" >&2
  echo "  <input.gz> will be replaced by <input>" >&2
  exit 11                ## distinctive non-0 exit code
}
input_gz="$1"            ## capture first argument "$1" as variable $input_gz

## print to stdout for posterity:

echo "input_gz:$input_gz"

## timestamp: to stdout e.g. '20220221114441' for 'Feb 21, 2022 at 11:44:41 am':

echo "begin:$(date +'%Y%m%d%H%M%S')"

## the actual work (capture stderr (2>&1) into stdout; and stdout into $result):

gunzip "$input_gz"

echo "finish:$(date +'%Y%m%d%H%M%S')"

