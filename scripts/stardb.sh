#!/usr/bin/env bash
#SBATCH -J stardb         ## $SLURM_JOB_NAME
#SBATCH -p compute        ## which set of physical compute servers to use
#SBATCH -q batch          ## quality of service (set of resource limits)
#SBATCH -t 72:00:00       ## time limit (max is 72h for '-q batch')
#SBATCH -c 16             ## numbers of CPUs (actually threads/vCPUs)
#SBATCH --mem 64g         ## RAM; 64 gigabyte
#SBATCH -o %x.%j.out      ## redirect stdout to $SLURM_JOB_NAME.$SLURM_JOB_ID.out
#SBATCH -e %x.%j.err      ## redirect stderr to $SLURM_JOB_NAME.$SLURM_JOB_ID.err

## settings:

cpus=16                   ## must match #SBATCH -c
index_n_bases=14          ## --genomeSAindexNbases [14]; set to min(14, log2(GenomeLength)/2 - 1)
container='/projects/researchit/crf/containers/crf_rnaseq.sif'

## check if exactly 2 arguments ($# has number of args); else print 
##   usage to stderr (>&2) and exit non-zero (error):

[ $# -eq 2 ] || {
  echo "Usage: sbatch stardb.sh <genome.fasta> <dir_out>" >&2
  exit 11
}
genome_fasta="$1"     ## load value of first argument "$1" into variable $genome_fasta
dir_out="$2"          ## load value of second argument "$2" into variable $dir_out

## print to stdout for posterity:

echo "genome_fasta:$genome_fasta"
echo "dir_out:$dir_out"
echo "index_n_bases:$index_n_bases"
echo "cpus:$cpus"
echo "container:$container"

module load singularity       ## make singularity command available

## star version captured in $version then printed to stdout:

version=$(singularity exec "$container" STAR --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:star_version: $version" >&2
  exit 21
}
echo "star_version:$version"

## timestamp: to stdout e.g. '20220221114441' for 'Feb 21, 2022 at 11:44:41 am':

echo "begin:$(date +'%Y%m%d%H%M%S')"

## the work (capture stderr (2>&1) into stdout; and stdout into $result):

result=$(singularity exec "$container" STAR \
  --runThreadN $cpus \
  --runMode genomeGenerate \
  --genomeDir "$dir_out" \
  --genomeSAindexNbases $index_n_bases \
  --genomeFastaFiles "$genome_fasta" 2>&1)

## $? is return value of last command; 0 if no error; if not 0,
##   signals error, then $result has error message from command:

[ $? -eq 0 ] || {
  echo "ERROR:star: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "star_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"

