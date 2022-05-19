#!/usr/bin/env bash
#SBATCH -J htseq_counts   ## $SLURM_JOB_NAME
#SBATCH -p compute        ## which set of physical compute servers to use
#SBATCH -q batch          ## quality of service (set of resource limits)
#SBATCH -t 72:00:00       ## time limit (max is 72h for '-q batch')
#SBATCH -c 24             ## numbers of CPUs (actually threads/vCPUs)
#SBATCH --mem 32g         ## RAM; g: gigabytes
#SBATCH -o %x.%j.out      ## redirect stdout to $SLURM_JOB_NAME.$SLURM_JOB_ID.out
#SBATCH -e %x.%j.err      ## redirect stderr to $SLURM_JOB_NAME.$SLURM_JOB_ID.err

## settings:

ncpus=24
mapq_min=10
mode='union'
nonunique='none'
secondary='ignore'
supplemental='ignore'
gtf='/projects/researchit/user15/db/GRCh38/release-103/gtf/Homo_sapiens.GRCh38.103.gtf'
out_file='htseq_counts.txt'
container='/projects/researchit/crf/containers/crf_rnaseq.sif'

## print to stdout for posterity:

echo "begin:$(date +'%Y%m%d%H%M%S')"
echo "ncpus:$ncpus"
echo "mapq_min:$mapq_min"
echo "mode:$mode"
echo "nonunique:$nonunique"
echo "secondary:$secondary"
echo "supplemental:$supplemental"
echo "gtf:$gtf"
echo "out_file:$out_file"
echo "container:$container"

module load singularity     ## make singularity available

## input alignment files:

files=($(ls *.bam))
echo "nfiles:${#files[@]}"
echo "files:${files[*]}"

## header row:

echo $(IFS=':' ; echo 'gene:'; echo "${files[*]}") | sed 's/Aligned.out.bam//g' | tr ':' $'\t' > "$out_file"

## the work (capture stderr (2>&1) into stdout; and stdout into $result):

result=$(singularity exec "$container" \
  htseq-count \
    -n $ncpus \
    -a $mapq_min \
    -m $mode \
    --nonunique=$nonunique \
    --secondary-alignments=$secondary \
    --supplementary-alignments=$supplemental \
    ${files[*]} "$gtf" >> "$out_file" 2>&1)

## $? is return value of last command; 0 if no error; if not 0,
##   signals error, then $result has error message from command:

[ $? -eq 0 ] || {
  echo "ERROR:htseq: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "htseq_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"

