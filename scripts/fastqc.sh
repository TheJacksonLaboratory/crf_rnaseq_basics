#!/usr/bin/env bash
#SBATCH -J fastqc        ## ## $SLURM_JOB_NAME
#SBATCH -p compute         ## partition
#SBATCH -q batch         ## qos
#SBATCH -c 2             ## number of threads must match 'cpus' variable
#SBATCH -t 72:00:00      ## time limit (72h)
#SBATCH --mem 8g         ## RAM (8 gigabytes)
#SBATCH -o %x.%j.out     ## stdout goes here ($SLURM_JOB_NAME.$SLURM_JOB_ID.out)
#SBATCH -e %x.%j.err     ## stderr goes here ($SLURM_JOB_NAME.$SLURM_JOB_ID.err)

## settings:

cpus=2                   ## must match #SBATCH -c
container=/projects/researchit/crf/containers/crf_rnaseq.sif

## make sure exactly one argument ($# has argument count)
##   or print usage to stderr (>&2) and exit non-zero (error):

[ $# -eq 1 ] || {
  echo "USAGE: sbatch fastqc.sh <reads.fastq>" >&2
  exit 11
}
fastq="$1"               ## set variable 'fastq' to value in first argument ($1)

echo "fastq:$fastq"
echo "container:$container"
echo "cpus:$cpus"

## make singularity command available:

module load singularity

## capture software version in $result then print to stdout:

result=$(singularity exec "$container" fastqc --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:fastqc_version: $result" >&2
  exit 21
}
echo "fastqc_version:$result"

## timestamp: prints out e.g. '20220221114441' for 'Feb 21, 2022 at 11:44:41 am'::

echo "begin: $(date +'%Y%m%d%H%M%S')"

## actual work (capture stderr (2>&1) into stdout; and stdout into $result):

result=$(singularity exec "$container" fastqc -t $cpus "$fastq" 2>&1)

## $? is return value of last command; 0 if no error; if not 0,
##   signals error, then $result has error message from command:

[ $? -eq 0 ] || {
  echo "ERROR:fastqc: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "fastqc_result:$result"

echo "finish: $(date +'%Y%m%d%H%M%S')"
