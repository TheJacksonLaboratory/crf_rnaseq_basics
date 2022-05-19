#!/usr/bin/env bash
#SBATCH -J feature_counts      ## $SLURM_JOB_NAME
#SBATCH -p compute             ## partition
#SBATCH -q batch               ## qos
#SBATCH -t 72:00:00            ## timelimit (72h)
#SBATCH -c 16                  ## number of cpus: must match 'cpus' variable below
#SBATCH --mem 16g              ## RAM (16 gigabytes)
#SBATCH -o %x.%j.out           ## stdout to $SLURM_JOB_NAME.$SLURM_JOB_ID.out
#SBATCH -e %x.%j.err           ## stderr to $SLURM_JOB_NAME.$SLURM_JOB_ID.err

## settings:

cpus=16                        ## must match SBATCH -c setting above
mapq_min=10                    ## minimum mapq score
min_overlap=5                  ## minimum exon overlap
container='/projects/researchit/crf/containers/crf_rnaseq.sif'

## check if exactly 2 arguments ($# has number of args); else print
##   usage to stderr (>&2) and exit non-zero (error):

[ $# -eq 2 ] || {
  echo "Usage: feature_counts.sh <reference.gtf> <output_file>" >&2
  exit 11
}
gtf="$1"                       ## load value of first argument "$1" into variable $gtf
output_file="$2"               ## load value of second argument "$2" into variable $output_file

## print to stdout for posterity:

echo "gtf:$gtf"
echo "output_file:$output_file"
echo "mapq_min:$mapq_min"
echo "min_overlap:$min_overlap"
echo "cpus:$cpus"
echo "container:$container"

## input alignment files:

files=($(ls *.bam))            ## outer '()' makes the returned value a list (splits on whitespace)
echo "nfiles:${#files[@]}"     ## prints the number of items in the list
echo "files:${files[*]}"       ## prints out the items in the list separated by spaces

module load singularity        ## make singularity command available

## timestamp: to stdout e.g. '20220221114441' for 'Feb 21, 2022 at 11:44:41 am':

echo "begin:$(date +'%Y%m%d%H%M%S')"

## the work (capture stderr (2>&1) into stdout; and stdout into $result):

result=$(singularity exec "$container" \
  featureCounts \
    -T $cpus \
    -Q $mapq_min \
    --minOverlap $min_overlap \
    --primary \
    --ignoreDup \
    -a "$gtf" \
    -o "$output_file" \
    ${files[*]} 2>&1)

## $? is return value of last command; 0 if no error; if not 0,
##   signals error, then $result has error message from command:

[ $? -eq 0 ] || {
  echo "ERROR:counts: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "counts_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
