#!/usr/bin/env bash
#SBATCH -J sras
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 1g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

sra_sh="$HOME/opt/crf_rnaseq_basics/sra.sh"
[ -e "$sra_sh" ] || {
  echo "ERROR: need to set 'sra_sh' to full path to 'sra.sh' file in repo directory 'crf_rnaseq_basics'." >&2
  exit 11
}

[ $# -eq 1 ] || {
  echo "USAGE: sbatch sras.sh <metadata.tsv>" >&2
  echo "  where <metadata.tsv> is tab delimited with sra_id in ... column." >&2
  echo "  creates .fastq.gz files corresponding to <sra_id>" >&2
  exit 13
}
metadata="$1"
echo "metadata:$metadata"

[ -f "$metadata" ] || {
  echo "metadata file '$metadata' is not a regular file." >&2
  exit 15
}

echo "begin:$(date +'%Y%m%d%H%M%S')"

base_url='http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo'

cut -d $'\t' -f 7 "$metadata" | while read url; do
  srx_id=${url##*=}
  srx_id=$(tr -d '"' <<< $srx_id)
  echo "srx_id:$srx_id"
  sra_ids=$(wget -qO- "$base_url&term=$srx_id" | grep $srx_id | cut -d ',' -f 1)
  for sra_id in "${sra_ids[@]}"; do
    echo -n "  ${sra_id}: "
    sbatch "$sra_sh" "$sra_id"
  done
done

echo "finish:$(date +'%Y%m%d%H%M%S')"

