#!/usr/bin/env bash
#SBATCH -J star           ## $SLURM_JOB_NAME
#SBATCH -p compute        ## which set of physical compute servers to use
#SBATCH -q batch          ## quality of service (set of resource limits)
#SBATCH -t 72:00:00       ## time limit (max is 72h for '-q batch')
#SBATCH -c 12             ## numbers of CPUs (actually threads/vCPUs)
#SBATCH --mem 64g         ## RAM; 64 gigabyte
#SBATCH -o %x.%j.out      ## redirect stdout to $SLURM_JOB_NAME.$SLURM_JOB_ID.out
#SBATCH -e %x.%j.err      ## redirect stderr to $SLURM_JOB_NAME.$SLURM_JOB_ID.err

## settings:

ncpus=12                   ## number of threads (should match SBATCH -c setting!!!)
genome_idx='star_idx'      ## star index directory
genome_gtf='reference.gtf'
container='/projects/researchit/crf/containers/crf_rnaseq.sif'
multimap_nmax=1            ## max number of genomic locations read maps to (else read unmapped)
mismatch_pmax=0.04         ## max proportion of read residues not matching template (else read unmapped)

## check if 2 or 3 arguments ($# has number of args); else print 
##   usage to stderr (>&2) and exit non-zero (error):

[ $# -eq 2 ] || [ $# -eq 3 ] || {
  echo "Usage: sbatch star.sh <prefix_out> <reads_forward.fastq> [<reads_reverse.fastq>]" >&2
  exit 11
}
prefix_out="$1"            ## output prefix (can include path); first argument is "$1"
reads_f="$2"               ## forward strand reads; second argument is "$3"
reads_r="$3"               ## optional reverse strand read for each forward read; 3d arg is "$4"

## print to stdout for posterity:

echo "prefix_out:$prefix_out"
echo "reads_f:$reads_f"
echo "reads_r:$reads_r"
echo "genome_idx:$genome_idx"
echo "genome_gtf:$genome_gtf"
echo "multimap_nmax:$multimap_nmax"
echo "mismatch_pmax:$mismatch_pmax"
echo "ncpus:$ncpus"
echo "container:$container"
echo "pwd:$(pwd)"

module load singularity       ## make singularity command available

## capture star version to $result, then print to stdout:

result=$(singularity exec "$container" STAR --version 2>&1)
[ $? -eq 0 ] || { 
  echo "ERROR:star_version: $result" >&2
  exit 21
}
echo "star_version:$result"

## timestamp: to stdout e.g. '20220221114441' for 'Feb 21, 2022 at 11:44:41 am':

echo "begin:$(date +'%Y%m%d%H%M%S')"

## the work (capture stderr (2>&1) into stdout; and stdout into $result):

result=$(singularity exec "$container" STAR \
  --runThreadN $ncpus \
  --genomeDir "$genome_idx" \
  --sjdbGTFfile "$genome_gtf" \
  --quantMode GeneCounts \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$prefix_out" \
  --readFilesIn $reads_f $reads_r >&2)

## $? is return value of last command; 0 if no error; if not 0,
##   signals error, then $result has error message from command:

[ $? -eq 0 ] || {
  echo "ERROR:star: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "star_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"

