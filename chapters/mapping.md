## Map reads

- [Home](../README.md)
- [Experimental Design](design.md)
- [Challenges to DE analysis](challenges.md)
- [Inputs](inputs.md)
- [Check sequence quality](fastqc.md)
- **[Map reads](mapping.md)** *(You are here)*
- [Generate expression matrix](count_matrix.md)
- [Import data into R](r_data.md)
- [Screen for outlier samples](outliers.md)
- [edgeR](edger.md)
- [DESeq2](deseq2.md)

---

### Table of contents

- [Introduction](#introduction)
- [Generate mapping index](#generate-mapping-index)
- [Read mapping](#read-mapping)
- [BAM and SAM formats](#bam-and-sam-formats)

---

### Introduction

The fundamental idea behind RNA-seq analysis is associating cDNA reads from a biological sample with identifiers of individual genes in order to estimate gene expression levels in that sample. This association is most accurately performed by using an intron-aware read aligner such as [STAR](https://github.com/alexdobin/STAR), [Bowtie2](https://github.com/BenLangmead/bowtie2), or [BWA-MEM](https://github.com/lh3/bwa). The job can be done with good accuracy much more quickly by using specialized programs developed for this purpose (sometimes called 'pseudo-aligners'), such as [kallisto](https://github.com/pachterlab/kallisto) and [Salmon](https://combine-lab.github.io/salmon/). These pseudo-aligners do not bother with the computationally intensive seed stitching step (described in the section below) full-fledged aligners perform. We demonstrate the use of the full-featured aligner STAR, because STAR can be used for a wide range of tasks, is very accurate, and fast enough to run in just a few minutes on our test dataset. We finish this chapter by describing the BAM and SAM sequence mapping formats.

[Table of contents](#table-of-contents)

---

### Generate mapping index

We will align the reads in our example dataset to the yeast genomic assembly using [STAR](https://github.com/alexdobin/STAR), whose documentation can be found [here](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). A description of the algorithm can be found [here](https://doi.org/10.1093/bioinformatics/bts635).

Briefly, the program works by building a [suffix array](https://en.wikipedia.org/wiki/Suffix_array) of the target set of sequences (in this case, contigs belonging to a genomic sequence assembly), sampling the sequences at a selectable interval (usually you should use the default settings). A suffix array is particularly efficient for searching for perfect matches of unpredictable length within a target string (in our case, segments of genomic sequence) in an amount of time proportional to the query (our short sequencing reads) length, rather than the longer target (genome sequence) length.

In the mapping phase of STAR, perfect matches starting at the 5' end of the read are used to search the suffix array. STAR also looks for perfect matches starting at adjustable intervals within the read. If part of the read contiguously matches the genomic sequence, but the rest does not, a new search is begun with the non-matching remainder of the read to see if there is a good match further downstream in the genomic sequence (which allows for interruptions in the alignment due to intervening introns). If the series of perfect matches (termed 'seeds') does not encompass the entire read, a 'stitching' phase is undertaken where a [dynamic programming](https://en.wikipedia.org/wiki/Dynamic_programming) algorithm similar to the relatively well-known [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) algorithm is used to find partially matching regions between each pair of seeds joining those seeds. The final stitched alignment is scored based on the number of matches versus the number of mismatches. The program is fairly insensitive to mismatches at the ends of the reads, in most cases obviating the need to trim adapter sequences from reads prior to mapping.

The basic command we will use for generating the mapping index (the set of suffix arrays and associated metadata) is:

```
STAR \
  --runThreadN $cpus \
  --runMode genomeGenerate \
  --genomeDir "$dir_out" \
  --genomeFastaFiles "$genome_fasta"
```

This command specifies the number of threads available `$cpus`, and the mode in which STAR is being run `genomeGenerate`, which triggers generation of a genome index. It also specifies the directory which will contain the index files `$dir_out`, and the input genomic sequence in FASTA format `$genome_fasta`.

In order to facilitate batch processing on the JAX HPC cluster, while automatically allocating a reasonable number of compute threads and RAM, we can use the included sbatch script, `stardb.sh`:

```
## make a separate directory in which to do the work for this step:

mkdir -p /fastscratch/$USER/rnaseq/stardb
cd /fastscratch/$USER/rnaseq/stardb

## make a link to the genomic assembly fasta file:

ln -s /fastscratch/$USER/rnaseq/reference/GCF_000146045.2_R64_genomic.fna .

## copy the stardb.sh script and edit it to set 'index_n_bases=10';
##   needed because yeast has a small genome; the defaults work better for
##   vertebrates; in general, it can be set to: min(14, log2(GenomeLength)/2 - 1):

cp ~/opt/crf_rnaseq_basics/scripts/stardb.sh .

sbatch stardb.sh GCF_000146045.2_R64_genomic.fna index_dir

## when the job completes, the output files are in the 'index_dir' directory specified:

$ ls -l index_dir
total 118528
-rw-r--r-- 1 user15 group12       122 Feb 11 10:14 chrLength.txt
-rw-r--r-- 1 user15 group12       327 Feb 11 10:14 chrNameLength.txt
-rw-r--r-- 1 user15 group12       205 Feb 11 10:14 chrName.txt
-rw-r--r-- 1 user15 group12       143 Feb 11 10:14 chrStart.txt
-rw-r--r-- 1 user15 group12  14680064 Feb 11 10:14 Genome
-rw-r--r-- 1 user15 group12       606 Feb 11 10:14 genomeParameters.txt
-rw-r--r-- 1 user15 group12      6180 Feb 11 10:14 Log.out
-rw-r--r-- 1 user15 group12 100296120 Feb 11 10:14 SA
-rw-r--r-- 1 user15 group12   6116787 Feb 11 10:14 SAindex
```

[Table of contents](#table-of-contents)

---

### Read mapping

We will map our reads to the index we just made using the STAR program again, but this time we will omit the `--runMode` setting, in which case the default value (`alignReads`) will be used, causing STAR to run in alignment mode. The basic command is:

```
STAR \
  --runThreadN $ncpus \
  --genomeDir "$genome_idx" \
  --sjdbGTFfile "$genome_gtf" \
  --quantMode GeneCounts \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$prefix_out" \
  --readFilesIn $reads_f $reads_r
```

Here we specify the number of threads `$ncpus`, the directory containing the STAR index `$genome_idx`, as well as the GTF2.2 formatted annotation file specifying the exon locations and their relationship to genes `$genome_gtf`. Setting `--quantMode GeneCounts` instructs SAM to generate an expression estimate for each gene. The format of the mapping file is specified with `--outSAMtype` to be a **BAM** (Binary Alignment Map) formatted file in which there is no sorting of mapping locations by coordinates (another option). We will discuss this format in more detail below. The resulting BAM file has nearby mappings from each read in a read pair (in the case of paired-end sequencing) guaranteed to be adjacent records in the BAM file. This condition is required by some downstream programs, such as `featureCounts`, which we use in the [expression matrix](count_matrix.md) chapter. We also specify the output filename prefix `$prefix_out`, the FASTQ file containing forward reads for one sample `$reads_f` and, if paired-end sequencing was performed, the FASTQ file containing the corresponding reverse reads `$reads_r`.

We can use the included convenience script `star.sh` to do the work on the JAX HPC cluster, but first we need to get our inputs together:

```
## make a directory for this step and go there:

mkdir -p /fastscratch/$USER/rnaseq/star
cd /fastscratch/$USER/rnaseq/star

## the script expects the STAR index to be in a subdirectory called 'star_idx'; we
##   can accomplish this by making a link to the `stardb.sh` output index directory:

ln -s /fastscratch/$USER/rnaseq/stardb/index_dir ./star_idx

## the script also expects the genomic annotation GTF2.2 file 'reference.gtf' in the
##   working directory; we can accomplish this by creating a link with the
##   expected name:

ln -s /projects/researchit/crf/data/reference/refseq_r64/GCF_000146045.2_R64_genomic.gtf reference.gtf

## finally, we link to the fastq sequences we want to map to the genome index:

ln -s /fastscratch/$USER/rnaseq/data/*.fastq .

## let's take a look:

$ ls -l
total 0
lrwxrwxrwx 1 user15 group12 93 Feb 11 10:15 reference.gtf -> /projects/researchit/crf/data/reference/refseq_r64/GCF_000146045.2_R64_genomic.gtf
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849386.fastq -> /fastscratch/user15/rnaseq/data/SRR11849386.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849387.fastq -> /fastscratch/user15/rnaseq/data/SRR11849387.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849388.fastq -> /fastscratch/user15/rnaseq/data/SRR11849388.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849389.fastq -> /fastscratch/user15/rnaseq/data/SRR11849389.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849390.fastq -> /fastscratch/user15/rnaseq/data/SRR11849390.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849391.fastq -> /fastscratch/user15/rnaseq/data/SRR11849391.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849392.fastq -> /fastscratch/user15/rnaseq/data/SRR11849392.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849393.fastq -> /fastscratch/user15/rnaseq/data/SRR11849393.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849394.fastq -> /fastscratch/user15/rnaseq/data/SRR11849394.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849395.fastq -> /fastscratch/user15/rnaseq/data/SRR11849395.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849396.fastq -> /fastscratch/user15/rnaseq/data/SRR11849396.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849397.fastq -> /fastscratch/user15/rnaseq/data/SRR11849397.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849398.fastq -> /fastscratch/user15/rnaseq/data/SRR11849398.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849399.fastq -> /fastscratch/user15/rnaseq/data/SRR11849399.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849400.fastq -> /fastscratch/user15/rnaseq/data/SRR11849400.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849401.fastq -> /fastscratch/user15/rnaseq/data/SRR11849401.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849402.fastq -> /fastscratch/user15/rnaseq/data/SRR11849402.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849403.fastq -> /fastscratch/user15/rnaseq/data/SRR11849403.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849404.fastq -> /fastscratch/user15/rnaseq/data/SRR11849404.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849405.fastq -> /fastscratch/user15/rnaseq/data/SRR11849405.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849406.fastq -> /fastscratch/user15/rnaseq/data/SRR11849406.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849407.fastq -> /fastscratch/user15/rnaseq/data/SRR11849407.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849408.fastq -> /fastscratch/user15/rnaseq/data/SRR11849408.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849409.fastq -> /fastscratch/user15/rnaseq/data/SRR11849409.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849410.fastq -> /fastscratch/user15/rnaseq/data/SRR11849410.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849411.fastq -> /fastscratch/user15/rnaseq/data/SRR11849411.fastq
lrwxrwxrwx 1 user15 group12 49 Feb 11 10:15 SRR11849412.fastq -> /fastscratch/user15/rnaseq/data/SRR11849412.fastq
lrwxrwxrwx 1 user15 group12 43 Feb 11 10:15 star_idx -> /fastscratch/user15/rnaseq/stardb/index_dir
```

Now that we have our inputs together, we can launch the script by sequentially feeding it one fastq file and a corresponding output prefix using the syntax:

```star.sh <prefix> <reads.fastq>```

If we had paired end data, we would use the syntax:

```star.sh <prefix> <forward.fastq> <reverse.fastq>```

We will use a slightly more complex bash script to manage the processing. Here we use the statement `p=${f%.fastq}`; the `${a%b}` notation indicates pattern substitution, and the `%` inside the curly braces means removal of the shortest match to the term on the right `b` within the variable specified on the left `$a`; so `${f%.fastq}` removes `.fastq` from the end of filename kept in the variable `$f`. We do this to provide a prefix for the outputs from `star.sh`. We do this so we can process all the same files within a single folder simultaneously without worrying about output files getting overwritten by different star jobs, and so we can link each set of outputs to the corresponding input FASTQ sequence read file:

```
## copy the script into the working directory:

cp ~/opt/crf_rnaseq_basics/scripts/star.sh .

## fire off the script once for each input file (demo data are not paired-end):

for f in *.fastq; do
  echo -n "$f: "         ## input fastq file
  p=${f%.fastq}          ## prefix: remove '.fastq' from end of file name
  sbatch star.sh $p $f
done
```

Once the batch jobs complete, we can get a quick idea of how well each job went by examining the corresponding `Log.final.out` file produced by STAR. The `uniquely mapped read number` and `uniquely mapped read %` are of particular interest:

```
$ cat SRR11849386Log.final.out
                                 Started job on |       Feb 11 10:21:44
                             Started mapping on |       Feb 11 10:21:45
                                    Finished on |       Feb 11 10:23:05
       Mapping speed, Million of reads per hour |       925.98

                          Number of input reads |       20577268
                      Average input read length |       74
                                    UNIQUE READS:
                   Uniquely mapped reads number |       17614016
                        Uniquely mapped reads % |       85.60%
                          Average mapped length |       73.84
                       Number of splices: Total |       143568
            Number of splices: Annotated (sjdb) |       130703
                       Number of splices: GT/AG |       139577
                       Number of splices: GC/AG |       368
                       Number of splices: AT/AC |       37
               Number of splices: Non-canonical |       3586
                      Mismatch rate per base, % |       0.15%
                         Deletion rate per base |       0.01%
                        Deletion average length |       1.35
                        Insertion rate per base |       0.00%
                       Insertion average length |       1.06
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       1393420
             % of reads mapped to multiple loci |       6.77%
        Number of reads mapped to too many loci |       253762
             % of reads mapped to too many loci |       1.23%
                                  UNMAPPED READS:
  Number of reads unmapped: too many mismatches |       0
       % of reads unmapped: too many mismatches |       0.00%
            Number of reads unmapped: too short |       1289618
                 % of reads unmapped: too short |       6.27%
                Number of reads unmapped: other |       26452
                     % of reads unmapped: other |       0.13%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%


```

[Table of contents](#table-of-contents)

---

### BAM and SAM formats

We requested that the alignment output from STAR be in BAM format, which is a binary encoding of the human-readable SAM format, documentated [here](https://samtools.github.io/hts-specs/SAMv1.pdf). The binary BAM format is good for saving space. However, we cannot see the alignment information easily. To get a better feel for what is in these alignment files (you don't normally need to do this, since the software we use in the next chapter accepts BAM formatted input), we can convert one of the output files into the human readable SAM format. To do so, we will use the [samtools](http://www.htslib.org/doc/samtools.html) program which is included in our singularity container. The conversion will output a new SAM file, while leaving the original BAM file unaffected:

```
## make samtools in a container available as the 'samtools' command:

module load singularity           ## load singularity module
alias samtools='singularity exec /projects/researchit/crf/containers/crf_rnaseq.sif samtools'

## try out our new alias:

$ samtools --version
samtools 1.11
Using htslib 1.11-4
Copyright (C) 2020 Genome Research Ltd.

## convert binary .bam format to human-readable .sam format;
##   -h keeps header; -o designates output file:

samtools view -h -o tmp.sam SRR11849401Aligned.out.bam           

## the top of the file has a header section, with each line starting w/ '@':

$ head tmp.sam
@HD     VN:1.4
@SQ     SN:NC_001133.9  LN:230218
@SQ     SN:NC_001134.8  LN:813184
@SQ     SN:NC_001135.5  LN:316620
@SQ     SN:NC_001136.10 LN:1531933
@SQ     SN:NC_001137.3  LN:576874
@SQ     SN:NC_001138.5  LN:270161
@SQ     SN:NC_001139.9  LN:1090940
@SQ     SN:NC_001140.6  LN:562643
@SQ     SN:NC_001141.2  LN:439888

## grep -v '^@': only pass lines that don't have '@' as first character;
##   that is, skip the header so we can see the alignment section:

$ grep -v '^@' tmp.sam | head
SRR11849401.496481      0       NC_001145.3     858406  255     74M     *       0       0       CTGGAAGCACCACCCACTTGCTGTGGCATACCGTTCGCATCGTAAGTCACAGCAGCACCAAACACAGGGAAATC    AAAAAEEEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEE      NH:i:1  HI:i:1AS:i:72 nM:i:0
SRR11849401.496482      0       NC_001144.5     454989  3       74M     *       0       0       GGAACATAGACAAGGAACGGCCCCAAAGTTGCCCTCTCCAAATTACAACTCGGGCACCGAAGGTACCAGATTTC    AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEE      NH:i:2  HI:i:1AS:i:72 nM:i:0
SRR11849401.496482      256     NC_001144.5     464126  3       74M     *       0       0       GGAACATAGACAAGGAACGGCCCCAAAGTTGCCCTCTCCAAATTACAACTCGGGCACCGAAGGTACCAGATTTC    AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEE      NH:i:2  HI:i:2AS:i:72 nM:i:0
SRR11849401.496484      0       NC_001144.5     416106  255     75M     *       0       0       GGCAACTTCCTTTTAAAGTGAATGATTCGCACAGCTCTGTCTTGTATAATTGTTTGGGAGTTTCCTGCACTTGGC   AA6AAEEE/EEEEEEEEE6EEEEEEE/AEEAEAEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<<     NH:i:1  HI:i:1AS:i:73 nM:i:0
SRR11849401.496485      0       NC_001137.3     259217  255     75M     *       0       0       CATAGTAGTCTGCTGTTTCACTTTAATAGCTTCAAATGGGCACAACATGATATCAGCGAGGAATTCAGCGGTCGC   AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEE     NH:i:1  HI:i:1AS:i:73 nM:i:0
SRR11849401.496489      16      NC_001137.3     267140  255     75M     *       0       0       ATAGCAGAAGCAGCACCAAGAATCATGGTGAAAAATAGAGGGAACGCTAGACCAGCAACTAGGGAGAAAAAAATC   EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA     NH:i:1  HI:i:1AS:i:73 nM:i:0
SRR11849401.496490      16      NC_001144.5     453440  3       66M     *       0       0       ATAAGACCCCATCTCCGGATAAACCAATTCCGGGGTGATAAGCTGTTAAGAAGAAAAGATAACTCC    EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA      NH:i:2  HI:i:1  AS:i:64 nM:i:0
SRR11849401.496490      272     NC_001144.5     462577  3       66M     *       0       0       ATAAGACCCCATCTCCGGATAAACCAATTCCGGGGTGATAAGCTGTTAAGAAGAAAAGATAACTCC    EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA      NH:i:2  HI:i:2  AS:i:64 nM:i:0
SRR11849401.496491      0       NC_001136.10    1491555 255     75M     *       0       0       ACCACCGTTAGACGATGAAGATGTAGACGTTGCCTTAGAAATTAGGCCTATGTTAGGCAAAGATGCAAAACATGT   AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE     NH:i:1  HI:i:1AS:i:73 nM:i:0
SRR11849401.496492      16      NC_001134.8     195480  255     75M     *       0       0       AAGAGGTCTATCGCCTAAGGAAAGAGCCCGTGAAATCATCAACAAGTGTGCTCATCCCGATTATCAAGCTTTGTT   EEEEEE<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA     NH:i:1  HI:i:1AS:i:73 nM:i:0
```

The top of the file is a header section, designated by marking the beginning of each line with a '`@`' character, and contains metadata relevant to the whole file. The alignments themselves come after the header and consist of at least (you can add optional fields) 11 tab-delimited fields. The first field is the sequence id of the query sequence (read) being aligned. The second field is an integer flag that encodes information about the alignment, such as whether it was a secondary or supplemental alignment. A secondary alignment is an alternative alignment of a largely overlapping portion of the same read to a different location in the genome. A supplemental alignment is an alternative alignment of the read to a different portion of the genome involving a largely non-overlapping portion of the same read. So the same fragment of a read aligning to different genomic locations results in secondary alignment, while alignment of non-overlapping fragments of the same read to different genomic locations results in supplemental alignments. See section 1.4 of the [SAM documentation](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information on the flag field. The third field is the sequence id of the target sequence to which the read was aligned. The fourth field is the position (1-based; that is counting begins with 1, not 0) of the left end of the mapping. The fifth field encodes the [MAPQ](http://maq.sourceforge.net/qual.shtml) or mapping quality score, which we will discuss in more depth below. The other columns are of less immediate interest for our workflow, but can be explored by perusing the documentation.

The **MAPQ score**, described in [this](https://genome.cshlp.org/content/18/11/1851.full.pdf) publication, is similar to the Phred score, in that it reflects the probability `P` that an error has been made. However, instead of referring to base-calling errors, the MAPQ score represents probabilities of incorrect mapping location. The way `P` is estimated differs for different types of mapping software, and is typically not well documented, so should be taken with a grain of salt. In any case, it should not be relied upon for comparing outputs from different programs. The MAPQ score usually reflects the number and (ideally) relative score of other alternative alignments of the same read, so if there are many other potential alignments of the read returned (to other regions of the genome), the MAPQ value of each alignment should reflect a high probability `P` of mapping error. Conversely, a uniquely mapping read will have a MAPQ value reflecting a very low probability of error. The formula for MAPQ score is similar to that used for Phred scores:

```
## R syntax:

## formula for Phred score (as a reminder);
##   P: probability of wrong base call:

Q = -10 * log10(P)

## formula for MAPQ;
##   P: probability of wrong alignment location:

MAPQ = round(-10 * log10(P))
```

For any read which maps to two positions, the lower scoring alignment position should have an estimated incorrect mapping probabilty `P` of no lower than `1 / 2`, which corresponds to a MAPQ value of no higher than:

```
> round(-10 * log10(1 / 2))
[1] 3
```

In general, we should be able to select for primary mappings by using a MAPQ cutoff of 3. In practice, this still implies a rather high probability of incorrect mapping (just under 50%), potentially introducing a lot of noise into our process (and why we don't just use the SAM/BAM flag field to screen out secondary and supplemental mappings). Typically, to get more reliable mappings of reads to genes, we like to include reads with an incorrect mapping probability no greater than 10% (or `P = 0.1`; which is still a pretty liberal cutoff), which would correspond to a lower MAPQ cutoff of `10`:

```
> round(-10 * log10(0.10))
[1] 10
```

We will use this cutoff in the following [expression matrix](count_matrix.md) chapter.

[Table of contents](#table-of-contents)

---

Next: [Generate expression matrix](count_matrix.md)  
Previous: [Check sequence quality](fastqc.md)  
