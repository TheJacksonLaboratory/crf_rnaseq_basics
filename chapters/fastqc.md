## Check sequence quality

- [Home](../README.md)
- [Experimental Design](design.md)
- [Challenges to DE analysis](challenges.md)
- [Inputs](inputs.md)
- **[Check sequence quality](fastqc.md)** *(You are here)*
- [Map reads](mapping.md)
- [Generate expression matrix](count_matrix.md)
- [Import data into R](r_data.md)
- [Screen for outlier samples](outliers.md)
- [edgeR](edger.md)
- [DESeq2](deseq2.md)

---

We will use the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) program to check for commonly encountered sequence quality issues. This program generates a long series of plots and tables summarizing various aspects of sequence quality for a FASTQ file. Some nice examples of output generated using good datasets as well as datasets with different types of issues are included on the [site](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and you are highly encouraged to go through those materials to get a better sense of how to interpret the output. More detailed descriptions of what each plot represents is found in the web pages listed [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

The basic command we will use for executing fastqc on each input FASTQ file will look like this:

```
fastqc -t $cpus "$fastq"
```

where `$cpus` is the number of vCPUS/threads available and `$fastq` is the path to the FASTQ file we want to evaluate.

The FastQC program and its dependencies are not available on the host system, like the `gunzip` command we used in the [Inputs](inputs.md) chapter. Instead, FastQC is packaged in the singularity container `/projects/researchit/crf/containers/crf_rnaseq.sif`, so the more complete command will look like this:

```
container=/projects/researchit/crf/containers/crf_rnaseq.sif
singularity exec "$container" fastqc -t $cpus "$fastq"
```

In order to facilitate non-interactive, multi-threaded batch execution on the JAX HPC cluster, using reasonable numbers of computational threads and RAM, we have included a slurm convenience script `fastqc.sh`, which wraps this process and can be invoked as follows:

```
sbatch fastqc.sh "$fastq"
```

If we have a directory of files, we can use this script and a bit of bash to process all the FASTQ files in the directory:

```
## make a separate directory in which to do the work for this step:

mkdir -p /fastscratch/$USER/rnaseq/fastqc

## copy the script over and link the data files:
cd /fastscratch/$USER/rnaseq
cp ~/opt/crf_rnaseq_basics/scripts/fastqc.sh fastqc  ## configure path for your install
ln -s $PWD/data/*.fastq fastqc
cd fastqc

$ ls -l
total 0
-rw-r--r-- 1 user15 group12 881 Feb 11 10:00 fastqc.sh
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849386.fastq -> /fastscratch/user15/rnaseq/data/SRR11849386.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849387.fastq -> /fastscratch/user15/rnaseq/data/SRR11849387.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849388.fastq -> /fastscratch/user15/rnaseq/data/SRR11849388.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849389.fastq -> /fastscratch/user15/rnaseq/data/SRR11849389.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849390.fastq -> /fastscratch/user15/rnaseq/data/SRR11849390.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849391.fastq -> /fastscratch/user15/rnaseq/data/SRR11849391.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849392.fastq -> /fastscratch/user15/rnaseq/data/SRR11849392.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849393.fastq -> /fastscratch/user15/rnaseq/data/SRR11849393.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849394.fastq -> /fastscratch/user15/rnaseq/data/SRR11849394.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849395.fastq -> /fastscratch/user15/rnaseq/data/SRR11849395.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849396.fastq -> /fastscratch/user15/rnaseq/data/SRR11849396.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849397.fastq -> /fastscratch/user15/rnaseq/data/SRR11849397.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849398.fastq -> /fastscratch/user15/rnaseq/data/SRR11849398.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849399.fastq -> /fastscratch/user15/rnaseq/data/SRR11849399.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849400.fastq -> /fastscratch/user15/rnaseq/data/SRR11849400.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849401.fastq -> /fastscratch/user15/rnaseq/data/SRR11849401.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849402.fastq -> /fastscratch/user15/rnaseq/data/SRR11849402.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849403.fastq -> /fastscratch/user15/rnaseq/data/SRR11849403.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849404.fastq -> /fastscratch/user15/rnaseq/data/SRR11849404.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849405.fastq -> /fastscratch/user15/rnaseq/data/SRR11849405.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849406.fastq -> /fastscratch/user15/rnaseq/data/SRR11849406.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849407.fastq -> /fastscratch/user15/rnaseq/data/SRR11849407.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849408.fastq -> /fastscratch/user15/rnaseq/data/SRR11849408.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849409.fastq -> /fastscratch/user15/rnaseq/data/SRR11849409.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849410.fastq -> /fastscratch/user15/rnaseq/data/SRR11849410.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849411.fastq -> /fastscratch/user15/rnaseq/data/SRR11849411.fastq
lrwxrwxrwx 1 user15 group12  49 Feb 11 10:00 SRR11849412.fastq -> /fastscratch/user15/rnaseq/data/SRR11849412.fastq

## launch jobs in parallel in the background; print corresponding job_ids:

for f in *.fastq; do
  echo -n "$f: "
  sbatch fastqc.sh "$f"
done
```

The main output file of interest will have the same prefix as the FASTQ file (everything up to '`.fastq`') followed by the extension '`_fastqc.html`'. This file can be viewed with your web browser and contains all the expected plots.

---

Next: [Map reads](mapping.md)  
Previous: [Inputs](inputs.md)  
