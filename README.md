# RNA-seq basics

## How to run a basic bulk RNA-seq analysis on the JAX HPC cluster

---

## Table of contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Requirements](#requirements)
- [Experimental data](#experimental-data)
- [Reference data](#reference-data)

## Links to subsequent chapters

- **[Home](README.md)** *(You are here)*
- [Experimental Design](./chapters/design.md)
- [Challenges to DE analysis](./chapters/challenges.md)
- [Inputs](./chapters/inputs.md)
- [Check sequence quality](./chapters/fastqc.md)
- [Map reads](./chapters/mapping.md)
- [Generate expression matrix](./chapters/count_matrix.md)
- [Import data into R](./chapters/r_data.md)
- [Screen for outlier samples](./chapters/outliers.md)
- [edgeR](./chapters/edger.md)
- [DESeq2](./chapters/deseq2.md)

---

## Introduction

This primer describes the steps involved in a basic bulk RNA-seq analysis. Here, bulk RNA-seq refers to a process in which RNA is extracted from a multi-cellular biological sample (such as a tissue, organ, whole organism, or collection of organisms), reverse transcribed into cDNA, and subjected to short-read sequencing (sequencing reads are shorter than 1000 bases in length). This process is most commonly used to sequence RNA present in each of several biological replicates from each of two or more experimental conditions of interest (e.g. treated vs. control; or wild-type vs mutant; male vs. female; etc.). These data are then analyzed to identify genes whose expression statistically (and hopefully substantively) differs between the conditions.

The data analysis involves mapping the reads to an annotated genome (other approaches may map directly to transcript sequences), then counting the number of reads that mapped to each annotated gene in each sample. This results in the generation of a rectangular matrix/table (called an **expression matrix**) with rows representing genes, columns representing samples, and each cell containing the number of reads from the column's sample that mapped to the row's gene. This part of the analysis is typically conducted in an HPC environment and requires several hours to compute. Once the expression matrix is in hand, we can analyze it for differential expression using conventional modeling approaches (t-test, anova, linear model) using specialized software that properly accounts for biological variability, the large number of genes being tested, and the relatively small number of biological replicates typically available. This differential expression analysis can usually be conducted on a modest laptop (at least 4 cores and 8 GB of RAM), but can also be carried out in the HPC environment.

[Table of contents](#table-of-contents)

---

## Overview

The procedural steps we cover in this tutorial are:

1. Setting up the input data for analysis.
2. Checking the quality of read sequences in fastq files using FASTQC.
3. Generating a STAR mapping index from a genomic sequence fasta file.
4. Mapping read sequences in FASTQ format to the STAR genome index, producing .bam output files.
5. Generating a raw expression matrix from the STAR mapping output.
6. Screening for outlier samples.
7. Differential expression analysis using edgeR.
8. Differential expression analysis using DESeq2.
  
Steps 1 through 5 are performed using bash commands. Subsequent steps are performed using R. Steps 1 through 5 are only presented in one form. By contrast, steps 7 and 8 represent different ways of doing essesially the same thing: testing for differential expression (DE). As a practical matter, our default choice for DE analysis is `DESeq2`, but the results from both methods tend to be similar, especially when replication is good and the size of effects associated with group membership are relatively strong. 

[Table of contents](#table-of-contents)

---

## Requirements

In order to understand this tutorial, you should have a **basic** understanding of the following subjects:

1. common bash commands
2. writing short bash scripts
3. common R commands
4. requesting resources from the SLURM job scheduler on the JAX HPC cluster
5. running a program installed in a singularity container.
6. familiarity with Welch's or Student's t-tests.
7. familiarity with the concept of a statistical distribution.

In addition, in order to understand the three-group and multivariate analyses sections in the `edgeR` and `DESeq2` chapters, familiarity with basic linear regression would be very helpful.

The scripts in this repository for steps 1 through 5 (up to and including expression matrix generation) assume you have at least 16 threads (or virtual CPUs), 64 GiB of RAM, and 200 GB of available storage on a compute node running some flavor of Linux. Steps 6 through 8 can be carried out on the JAX HPC cluster but can (unless your data are really big) more conveniently carried out on your laptop.

For convenience, we have included the the expression matrix `feature_counts.txt` and the metadata file `metadata.GSE151185.txt` in the `data` subdirectory of this repo. You can use these two files to skip steps 1 through 5, and proceed directly to steps 6 through 8.

On the JAX HPC cluster, you can download this repository (including the scripts) into a directory (say, for instance `opt`) under your `/home` directory by changing to the parent directory and entering:

```
mkdir -p ~/opt
cd ~/opt
git clone https://github.com/TheJacksonLaboratory/crf_rnaseq_basics.git
```

This will create the subdirectory `~/opt/crf_rnaseq_basics`, containing the contents of this repository.

On the JAX HPC cluster, the dependencies for all the steps (including R) are present in the Singularity container:

```
/projects/researchit/crf/containers/crf_rnaseq.sif
```

A list of software versions in the container can be obtained by running:

```
module load singularity
singularity run /projects/researchit/crf/containers/crf_rnaseq.sif
```

If you decide to do steps 6 through 8 on your laptop, you will need the statistical computing software `R` along with the `Bioconductor` packages `edgeR`, `limma`, and `DESeq2`. Once you have R installed on your laptop (may require elevated privileges), you can install the required R packages locally (should not require elevated privileges). To do so, after you launch the R console or RStudio, you can enter the following commands:

```
## First see if you have Bioconductor's package manager installed, and 
##   install it if you don't. The 'install.packages()' command might
##   queue you to select a package repository. You can pick anywhere 
##   geographically nearby. The 'require()' function returns FALSE if
##   the requested package cannot be loaded (because it is not
##   currently installed):

if (!require("BiocManager", quietly=T)) install.packages("BiocManager")

## Now we can use Bioconductor's package manager to install the
##   needed Bioconductor packages and their dependencies. The 
##   '::' indicates that the 'install()' function requested should
##   come from the BiocManager package:

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("DESeq2")
```

[Table of contents](#table-of-contents)

---

## Experimental data

The experiment we will be looking at was conducted in the model organism *Saccharomyces cerevisiae* (baker's/ale yeast). We picked this organism because the datasets and associated genome are smaller than typical for vertebrates, and therefore require fewer resources for demonstrating computational procedures.

These test data are derived from [GEO](https://www.ncbi.nlm.nih.gov/geo/) dataset [GSE151185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151185). We have made the sequence read files available on the JAX HPC cluster in the `/projects/researchit/crf/data/rnaseq/GSE151185` directory. This folder also includes a metadata file we made called `metadata.GSE151185.txt` which describes each biological sample and connects it to the corresponding sequencing read file. We describe the metadata file in more detail, as well as how to copy and setup the data for analysis in the chapter on [Setting up inputs](chapters/inputs.md).

These data were generated for a [study](https://dx.doi.org/10.1074%2Fjbc.RA120.015402) looking at the mechanisms by which caloric restriction increases longevity in yeast. Caloric restriction has been shown to extend the average lifespan of virtually every species of eukaryotic model organism reared in the laboratory (a situation in which calories tend to be much more readily available than the natural environment). The investigators noted that conditioned media (media collected from growing cultures by filtering out the cells) from calorie restricted cultures at stationary phase was able to extend the 'chronological longevity' of fresh yeast cultures, even when supplemented to reach normal calorie levels. This suggested a factor released from calorie restricted yeast was able to increase longevity in the absence of caloric restriction. Chronological longevity refers to the amount of time that yeast continue to be able to proliferate after entering stationary phase. This is tested by taking samples from stationary phase cultures at different timepoints and seeing how well they grow when transferred to fresh media.

In order to investigate the mechanism of longevity extension by conditioned media from calorie restricted cultures, researchers exposed yeast to three types of media, or 'treatments': 

1. Normal media supplemented with conditioned media from calorie restricted cultures (`CR`).
2. Normal media supplemented with conditioned media from non-calorie restricted cultures (`NR`).
3. Normal media supplemented with water (`water`).

Each treatment was represented by three timepoints (`0` hours, `24` hours and `96` hours). Each timepoint:treatment combination was represented by three biological replicates (where a replicate is a separate culture container). Cells were collected at the appropriate timepoint, mRNA was isolated, cDNA was prepared and then sequenced. After differential expression analysis, the investigators highlighted the gene `YCL064C` as being a likely player in the mechanism by which treatment `CR` extended lifespan.

[Table of contents](#table-of-contents)

## Reference data

The main idea behind RNA-seq is to associate cDNA sequencing reads from a sample with the corresponding gene, then use the number of reads associated with different genes to estimate relative gene expression. In practice, this association of reads with genes is done by searching for sequence matches between each read and the genomic sequence. Then genomic annotation containing the coordinates of each gene's exons is cross-referenced against the locations in the genomic sequence matching the read. If the read's matching location overlaps a gene's exon location, the read is added to the read count tally for that gene. In order to carry out this process we need the genomic assembly sequence and associated exon annotations, which include exon coordinates and relationship to gene ids. 

The reference genome and associated annotation for well-known model-organisms (including human) will usually come from either [Gencode](https://www.gencodegenes.org/) or [Refseq](https://www.ncbi.nlm.nih.gov/refseq/). The two sources provide the same genomic assembly sequences, but the identifiers for these sequences can differ between sources. The main difference, is that Gencode annotation is a superset of Refseq annotation. Refseq only annotates features which they deem have been identified with high confidence. Gencode supplements these with lower-confidence gene models, for instance identified by *ab initio* gene prediction programs. 

For our example analysis, we found our genome by searching for `Saccharomyces cerevisiae` [here](https://www.ncbi.nlm.nih.gov/genome), which led us to the reference genome [here](https://www.ncbi.nlm.nih.gov/genome/15?genome_assembly_id=22535). Following the `Download sequence and annotation from RefSeq` link led to the FTP folder [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/), which contained the files we will need.

---

Next: [Experimental design](./chapters/design.md)
