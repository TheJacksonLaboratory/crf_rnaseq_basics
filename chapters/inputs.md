## Inputs

- [Home](../README.md)
- [Experimental Design](design.md)
- [Challenges to DE analysis](challenges.md)
- **[Inputs](inputs.md)** *(You are here)*
- [Check sequence quality](fastqc.md)
- [Map reads](mapping.md)
- [Generate expression matrix](count_matrix.md)
- [Import data into R](r_data.md)
- [Screen for outlier samples](outliers.md)
- [edgeR](edger.md)
- [DESeq2](deseq2.md)
---

### Table of contents

- [Introduction](#introduction)
- [Data setup](#data-setup)
- [File decompression](#file-decompression)
- [Metadata](#metadata)
- [FASTA format](#fasta-format)
- [GTF2 format](#gtf2-format)
- [FASTQ format](#fastq-format)

---

### Introduction

This chapter will describe the input files we use for the rest of the analysis. It walks through how to set up working directories, how to transfer the data, and how to decompress the data. Subsequent sections in this chapter go on to describe the format of each file, including the metadata file, the FASTA formatted genomic sequence, the GTF2 formatted genomic annotation, and finally the FASTQ formatted experimental sequencing reads.

In order to perform the example analysis, you will need to complete the steps in the first two sections, [Data setup](#data-setup) and [File decompression](#file-decompression). If you are familiar with the file formats, you can safely skip the subsequent sections.

[Table of contents](#table-of-contents)

---

### Data setup

The experiment we will be looking at was conducted in the model organism *Saccharomyces cerevisiae* (baker's/ale yeast). The usual starting point for an RNA-seq analysis is a set of **FASTQ** formatted files containing reads, a meta-data text file mapping file names to conditions (as well other variables of potential interest), a **FASTA** file containing genomic sequence, and a **GTF** file containing exon coordinates and associations with gene ids. The FASTQ files are typically created by the sequencing facility. Often, they may be compressed, in which case it may (depending on which software you are using) be necessary to decompress the files prior to further use. Metadata files are typically tab-delimited, with one column containing the FASTQ file name, and other columns indicating the associated condition as well as values for other variables of interest.

Our example input sequence read data and a corresponding metadata file can be found in:

```
/projects/researchit/crf/data/rnaseq/GSE151185
```

We will get started by copying the data over to a working directory we create for this purpose under `/fastscratch`:

```
echo $USER                                 ## take a look at the predefined USER variable
mkdir -p /fastscratch/$USER/rnaseq/data
cd /fastscratch/$USER/rnaseq/data
rsync -v /projects/researchit/crf/data/rnaseq/GSE151185/* .

## our fastq files are compressed (can tell from '.gz' suffix):

$ ls -l
total 24918400
-r--r--r-- 1 user15 group12       1445 Feb 11 09:27 metadata.GSE151185.txt
-r--r--r-- 1 user15 group12  915317492 Feb 11 09:25 SRR11849386.fastq.gz
-r--r--r-- 1 user15 group12 1019003063 Feb 11 09:25 SRR11849387.fastq.gz
-r--r--r-- 1 user15 group12 1795906600 Feb 11 09:25 SRR11849388.fastq.gz
-r--r--r-- 1 user15 group12  796989960 Feb 11 09:25 SRR11849389.fastq.gz
-r--r--r-- 1 user15 group12  669081323 Feb 11 09:26 SRR11849390.fastq.gz
-r--r--r-- 1 user15 group12 1086323700 Feb 11 09:26 SRR11849391.fastq.gz
-r--r--r-- 1 user15 group12 1285240997 Feb 11 09:26 SRR11849392.fastq.gz
-r--r--r-- 1 user15 group12  762600699 Feb 11 09:26 SRR11849393.fastq.gz
-r--r--r-- 1 user15 group12 1200210880 Feb 11 09:26 SRR11849394.fastq.gz
-r--r--r-- 1 user15 group12  678294178 Feb 11 09:26 SRR11849395.fastq.gz
-r--r--r-- 1 user15 group12  563999341 Feb 11 09:26 SRR11849396.fastq.gz
-r--r--r-- 1 user15 group12  440662866 Feb 11 09:26 SRR11849397.fastq.gz
-r--r--r-- 1 user15 group12  528612140 Feb 11 09:26 SRR11849398.fastq.gz
-r--r--r-- 1 user15 group12  400186401 Feb 11 09:26 SRR11849399.fastq.gz
-r--r--r-- 1 user15 group12 5881481616 Feb 11 09:26 SRR11849400.fastq.gz
-r--r--r-- 1 user15 group12  467012406 Feb 11 09:26 SRR11849401.fastq.gz
-r--r--r-- 1 user15 group12 1763551773 Feb 11 09:26 SRR11849402.fastq.gz
-r--r--r-- 1 user15 group12  673537185 Feb 11 09:26 SRR11849403.fastq.gz
-r--r--r-- 1 user15 group12  411422497 Feb 11 09:26 SRR11849404.fastq.gz
-r--r--r-- 1 user15 group12  716508863 Feb 11 09:26 SRR11849405.fastq.gz
-r--r--r-- 1 user15 group12  397152943 Feb 11 09:27 SRR11849406.fastq.gz
-r--r--r-- 1 user15 group12  593729749 Feb 11 09:27 SRR11849407.fastq.gz
-r--r--r-- 1 user15 group12  618690779 Feb 11 09:27 SRR11849408.fastq.gz
-r--r--r-- 1 user15 group12  479065857 Feb 11 09:27 SRR11849409.fastq.gz
-r--r--r-- 1 user15 group12  426470532 Feb 11 09:27 SRR11849410.fastq.gz
-r--r--r-- 1 user15 group12  406697846 Feb 11 09:27 SRR11849411.fastq.gz
-r--r--r-- 1 user15 group12  537321799 Feb 11 09:27 SRR11849412.fastq.gz
```

We copied the `fastq` files containing the sequencing reads from our experiment, so as not to lose or accidentally overwrite the original versions. We'll have to decompress these files, which will overwrite our compressed copies, but leave our originals unaffected.

We will avoid making actual copies of the reference genome, and instead make links to save space:

```
mkdir /fastscratch/$USER/rnaseq/reference
cd /fastscratch/$USER/rnaseq/reference

## make link to genomic sequences:

ln -s /projects/researchit/crf/data/reference/refseq_r64/GCF_000146045.2_R64_genomic.fna .

## make link to genomic annotation:

ln -s /projects/researchit/crf/data/reference/refseq_r64/GCF_000146045.2_R64_genomic.gtf .

$ ls -l
total 0
lrwxrwxrwx 1 user15 group12 93 Feb 11 09:31 GCF_000146045.2_R64_genomic.fna -> /projects/researchit/crf/data/reference/refseq_r64/GCF_000146045.2_R64_genomic.fna
lrwxrwxrwx 1 user15 group12 93 Feb 11 09:31 GCF_000146045.2_R64_genomic.gtf -> /projects/researchit/crf/data/reference/refseq_r64/GCF_000146045.2_R64_genomic.gtf
```

[Table of contents](#table-of-contents)

---

### File decompression

As we saw in the previous section, our FASTQ files, which contain read sequences from our example experiment are compressed. The `.gz` extension indicates they can be decompressed using the `gunzip` command. The basic syntax is `gunzip <input.gz>`. Executing this command will cause `<input.gz>` to be replaced by a file `<input>` (same name, but no `.gz` suffix), containing the decompressed data. This process can take a while. In order to help speed things along we will use a convenience slurm script `gunzip.sh` from this repo in order to parallelize the task. The structure of this script is very similar to the structure of other convenience scripts included in this repo, so we'll review it here. The `gunzip.sh` script has a header in it specifying the program to use when interpreting the script and the resources to be requested from the [Slurm](https://slurm.schedmd.com/) workload manager:

```
#!/usr/bin/env bash
#SBATCH -J gunzip         ## name of the job in slurm queueing system
#SBATCH -p compute        ## which set of physical compute servers to use
#SBATCH -q batch          ## quality of service (set of resource limits)
#SBATCH -t 72:00:00       ## time limit (max is 72h for '-q batch')
#SBATCH -c 1              ## numbers of CPUs (actually threads/vCPUs)
#SBATCH --mem 1g          ## RAM; g: gigabytes
#SBATCH -o '%x.%j.out'    ## redirect stdout to $SLURM_JOB_NAME.$SLURM_JOB_ID.out
#SBATCH -e '%x.%j.err'    ## redirect stderr to $SLURM_JOB_NAME.$SLURM_JOB_ID.err
```

Each line in the header begins with a `#` character, which causes `bash` to consider the rest of the line a comment. This means that the script can be executed with whatever resources your terminal session currently has available, without invoking Slurm. In general, the syntax for invoking such a script with and without Slurm is:

```
## [<arg1> [<arg2> ... ]] means a list of 0 or more arguments:

sbatch <script.sh> [<arg1> [<arg2> ... ]]         ## invoke with SLURM
bash <script.sh> [<arg1> [<arg2> ... ]]           ## invoke without SLURM
```

This header is followed by by a section which processes the command line arguments. Our script is set to look for one command line argument, and returns a usage message if that condition is not met. Otherwise it loads the first (only) argument into the variable `input_gz`. As a reminder: `$#` is the number of arguments, the first of which is stored in `$1`, the second in `$2`, etc. The `[ $# -eq 1 ]` tests for integer equality. The `||` is a logical 'or' operator, which means the part on the right is only executed if the part on the left is false. Normal execution returns an exit code of 0, while any other exit code indicates an error of some sort.

```
## $# (number of arguments) equals 1 or print usage to stderr (>&2) and exit non-0:

[ $# -eq 1 ] || {
  echo "Usage: sbatch gunzip.sh <input.gz>" >&2
  echo "  <input.gz> will be replaced by <input>" >&2
  exit 11                ## distinctive non-0 exit code
}
input_gz="$1"            ## capture first argument as variable $input_gz
```

The argument processing section is followed by a line recording the input file in the standard output capture file (one of our two log files) `echo "input_gz:$input_gz"`.

We then get to the actual decompression command, capturing the beginning and ending times with a command of the form `echo "begin:$(date +'%Y%m%d%H%M%S')"`. This will record the datetime in an all digit format with 4-digit year `%Y`, 2 digits each for month `%m`, day `%d`, hour `%H`, minute `%M`, and second `%S`. So `Feb 21, 2022 at 11:44:41 am` will be recorded as `20220221114441`:

```
echo "begin:$(date +'%Y%m%d%H%M%S')"
gunzip "$input_gz"
echo "finish:$(date +'%Y%m%d%H%M%S')"
```

Now that we understand the structure of our scripts, we'll go ahead and decompress our FASTQ files:

```
## change to our working directory:
cd /fastscratch/$USER/rnaseq/data

## a way to loop across files with the same file-extension (.gz):

for f in *.gz; do
  echo $f
done

## using that loop, let's see what our command looks like before we execute (so we don't mess up):

for f in *.gz; do                                      ## f set to successive files matching '*.gz'
  echo -n "$f:"                                        ## print filename without newline
  echo "sbatch ~/opt/crf_rnaseq_basics/scripts/gunzip.sh $f"   ## print the command we intend to execute
done

## if everything looks right, decompress all the compressed files:

for f in *.gz; do                                      ## f set to successive files matching '*.gz'
  echo -n "$f:"                                        ## print filename without newline
  sbatch ~/opt/crf_rnaseq_basics/scripts/gunzip.sh "$f"        ## gunzip file; prints slurm jobid (w/ newline)
done

## we can periodically check our progress:

squeue -u $USER

## or monitor it continuously in a dedicated terminal window, refeshing every 5 seconds:

while [ 0 ]; do
  squeue -u $USER
  echo
  sleep 5
done

## or equivalently:

while [ 0 ]; do squeue -u $USER; echo; sleep 5; done

## when everything finishes (jobs no longer listed), we'll check error file size (should be zero):

$ ls -l *.err
-rw-r--r-- 1 user15 group12 0 Feb 11 09:48 gunzip.13443241.err
-rw-r--r-- 1 user15 group12 0 Feb 11 09:48 gunzip.13443242.err
-rw-r--r-- 1 user15 group12 0 Feb 11 09:48 gunzip.13443243.err

...   ## skip a bunch of 0 size files
-rw-r--r-- 1 user15 group12 0 Feb 11 09:48 gunzip.13443266.err
-rw-r--r-- 1 user15 group12 0 Feb 11 09:48 gunzip.13443267.err

## and make sure everything finished with an exit status of zero;
##   grep returns lines matching pattern in first argument; '^' is
##   beginning of line, so '^State:' only matches 'State:' at the beginning
##   of a line:

$ grep '^State:' *.out
gunzip.13443241.out:State: COMPLETED (exit code 0)
gunzip.13443242.out:State: COMPLETED (exit code 0)
gunzip.13443243.out:State: COMPLETED (exit code 0)
...   ## skip a bunch of 0 exit codes
gunzip.13443266.out:State: COMPLETED (exit code 0)
gunzip.13443267.out:State: COMPLETED (exit code 0)

## clean up logs:

mkdir logs
mv *.out *.err logs

## the compressed versions are replaced by the uncompressed versions,
##   which are typically 2-4x larger in size:

$ ls -l
total 129537152
drwxr-xr-x 2 user15 group12        4096 Feb 11 09:53 logs
-r--r--r-- 1 user15 group12        1445 Feb 11 09:27 metadata.GSE151185.txt
-r--r--r-- 1 user15 group12  4729941466 Feb 11 09:25 SRR11849386.fastq
-r--r--r-- 1 user15 group12  5249931756 Feb 11 09:25 SRR11849387.fastq
-r--r--r-- 1 user15 group12 10172674994 Feb 11 09:25 SRR11849388.fastq
-r--r--r-- 1 user15 group12  4036599550 Feb 11 09:25 SRR11849389.fastq
-r--r--r-- 1 user15 group12  3409229486 Feb 11 09:26 SRR11849390.fastq
-r--r--r-- 1 user15 group12  5599969862 Feb 11 09:26 SRR11849391.fastq
-r--r--r-- 1 user15 group12  6810613644 Feb 11 09:26 SRR11849392.fastq
-r--r--r-- 1 user15 group12  4204589032 Feb 11 09:26 SRR11849393.fastq
-r--r--r-- 1 user15 group12  6735960046 Feb 11 09:26 SRR11849394.fastq
-r--r--r-- 1 user15 group12  3411335794 Feb 11 09:26 SRR11849395.fastq
-r--r--r-- 1 user15 group12  2843857938 Feb 11 09:26 SRR11849396.fastq
-r--r--r-- 1 user15 group12  2365506634 Feb 11 09:26 SRR11849397.fastq
-r--r--r-- 1 user15 group12  2634588062 Feb 11 09:26 SRR11849398.fastq
-r--r--r-- 1 user15 group12  1995974090 Feb 11 09:26 SRR11849399.fastq
-r--r--r-- 1 user15 group12 29775508460 Feb 11 09:26 SRR11849400.fastq
-r--r--r-- 1 user15 group12  2724690290 Feb 11 09:26 SRR11849401.fastq
-r--r--r-- 1 user15 group12  9034970912 Feb 11 09:26 SRR11849402.fastq
-r--r--r-- 1 user15 group12  3466868914 Feb 11 09:26 SRR11849403.fastq
-r--r--r-- 1 user15 group12  2070082360 Feb 11 09:26 SRR11849404.fastq
-r--r--r-- 1 user15 group12  3650270142 Feb 11 09:26 SRR11849405.fastq
-r--r--r-- 1 user15 group12  2223135858 Feb 11 09:27 SRR11849406.fastq
-r--r--r-- 1 user15 group12  3013827614 Feb 11 09:27 SRR11849407.fastq
-r--r--r-- 1 user15 group12  3124481144 Feb 11 09:27 SRR11849408.fastq
-r--r--r-- 1 user15 group12  2396916428 Feb 11 09:27 SRR11849409.fastq
-r--r--r-- 1 user15 group12  2189989816 Feb 11 09:27 SRR11849410.fastq
-r--r--r-- 1 user15 group12  2051538590 Feb 11 09:27 SRR11849411.fastq
-r--r--r-- 1 user15 group12  2721115192 Feb 11 09:27 SRR11849412.fastq
```

[Table of contents](#table-of-contents)

---

### Metadata

We can take a look at the top of the metadata file to get a better feel for its format, as well as what this dataset is about:

```
cd /fastscratch/$USER/rnaseq/data

## take a look at first 10 (head's default) lines:

$ head metadata.GSE151185.txt
Sample        SRX           SRR           Label        treatment  time   rep
GSM4568115    SRX8399630    SRR11849386   Waterlog_1   water      0      1
GSM4568116    SRX8399631    SRR11849387   Waterlog_2   water      0      2
GSM4568117    SRX8399632    SRR11849388   Waterlog_2   water      0      3
GSM4568118    SRX8399606    SRR11849389   Water24_1    water      24     1
GSM4568119    SRX8399607    SRR11849390   Water24_2    water      24     2
GSM4568120    SRX8399608    SRR11849391   Water24_3    water      24     3
GSM4568121    SRX8399609    SRR11849392   Water96_1    water      96     1
GSM4568122    SRX8399610    SRR11849393   Water96_2    water      96     2
GSM4568123    SRX8399611    SRR11849394   Water96_3    water      96     3
```

In the analysis shown below, the main columns/variables of interest will be the `treatment` and `time`.  The `SRR` column can easily be converted into a filename by adding the suffix `.fastq`. The `Sample` and `SRX` columns contain other identifiers from the Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) and Sequence Read Archive ([SRA](https://www.ncbi.nlm.nih.gov/sra)) data repositories at [NCBI](https://www.ncbi.nlm.nih.gov/home/about/). The data we are working with correspond to [this](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151185) experiment.

[Table of contents](#table-of-contents)

---

### FASTA format

The reference genome will be represented by a sequence file in FASTA format, and an annotation file in GTF format. We'll look at both formats below. The sequences themselves are generally the same between comparable versions of Gencode and Refseq, but the sequence identifiers can differ. The annotation files from Gencode will typically specify many more genes and trascripts than Refseq, reflecting the more conservative approach to gene annotation used by the Refseq effort, which does not include purely hypothetical (predicted) genes included in Gencode.

The [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format is a relatively simple format for sequence data. We can take a look at one by entering:

```
cd /fastscratch/$USER/rnaseq

$ head reference/GCF_000146045.2_R64_genomic.fna
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
ccacaccacacccacacacccacacaccacaccacacaccacaccacacccacacacacacatCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
TCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATC
CAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATACTGTTCTTCTACCCACCATAT
TGAAACGCTAACAAATGATCGTAAATAACACACACGTGCTTACCCTACCACTTTATACCACCACCACATGCCATACTCAC
CCTCACTTGTATACTGATTTTACGTACGCACACGGATGCTACAGTATATACCATCTCAAACTTACCCTACTCTCAGATTC
CACTTCACTCCATGGCCCATCTCTCACTGAATCAGTACCAAATGCACTCACATCATTATGCACGGCACTTGCCTCAGCGG
TCTATACCCTGTGCCATTTACCCATAACGCCCATCATTATCCACATTTTGATATCTATATCTCATTCGGCGGTcccaaat
attgtataaCTGCCCTTAATACATACGTTATACCACTTTTGCACCATATACTTACCACTCCATTTATATACACTTATGTC
```

The first character on the first line is a `>` character, which indicates that line is a 'description line' containing a sequence identifier and description for the sequence. The first item (delimited by whitespace) to appear on a description line is the required sequence identifier. In this case the first sequence identifier is `NC_001133.9`, which corresponds to chromosome 1 of the Baker's Yeast genome. Any additional words are optional, and constitute the description of the sequence. Lines occuring after a description line are the nucleotide sequence itself, which can optionally be wrapped across multiple lines. Each sequence record is delimited by the next description line, or by the end of the file. The use of lowercase indicates the sequence has been 'soft-masked' which is sometimes done with repetitive sequence in order to prevent spurious sequence alignments. Many sequence aligners will not initiate sequence alignments in such a region, but will extend alignments started in adjacent sequence through such regions. Without this feature, the sequence aligner can initiate an enormous number of partial alignments based on fragmentary matches of some reads to repetitive sequence scattered throughout the genome.

[Table of contents](#table-of-contents)

---

### GTF2 format

The reference genome annotation needed is in the [GTF2.2](https://mblab.wustl.edu/GTF22.html) format. This is a tab-delimited file with 9 columns. The first few lines are preceded by `#` comment characters and constitute a header section, describing the format version and other meta-data. Each line after the header section represents a single genomic feature. The first column is a sequence identifier that matches a sequence identifier in the corresponding reference genome FASTA sequence file. The second column is the source of the annotation line. The third column represents the type of feature. The most important type for bulk RNA-seq analysis is `exon` records. The fourth and fifth columns represent the `start` and `end` positions of the feature (one based, closed; that is counting begins at 1; `start` and `end` positions are included in the feature, so a feature that encompasses a single base will have `start == end`). The sixth column is the score (empty/undetermined values in any column are designated with a `.` symbol), the seventh is the strand on which the feature occurs (one of `+`, `-`, or `.`). The eighth column is the open reading frame (only matters for protein features), and the ninth is a collection of attributes in `key1 "value1"; key2 "value2";` format. The most important attribute for our analysis is the `gene_id` associated with each `exon` record. In RNA-seq analysis we generate expression estimates for each gene by mapping reads to the genomic sequence, then counting the number of reads whose mappings overlap with exons associated with that gene:

```
$ head reference/GCF_000146045.2_R64_genomic.gtf
#gtf-version 2.2
#!genome-build R64
#!genome-build-accession NCBI_Assembly:GCF_000146045.2
#!annotation-source SGD R64-3-1
NC_001133.9 RefSeq  gene    1807    2169    .   -   .   gene_id "YAL068C"; transcript_id ""; db_xref "GeneID:851229"; gbkey "Gene"; gene "PAU8"; gene_biotype "protein_coding"; locus_tag "YAL068C"; partial "true";
NC_001133.9 RefSeq  transcript  1807    2169    .   -   .   gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "GeneID:851229"; gbkey "mRNA"; gene "PAU8"; locus_tag "YAL068C"; partial "true"; product "seripauperin PAU8"; transcript_biotype "mRNA";
NC_001133.9 RefSeq  exon    1807    2169    .   -   .   gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "GeneID:851229"; gene "PAU8"; locus_tag "YAL068C"; partial "true"; product "seripauperin PAU8"; transcript_biotype "mRNA"; exon_number "1";
NC_001133.9 RefSeq  CDS 1810    2169    .   -   0   gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "SGD:S000002142"; db_xref "GeneID:851229"; experiment "EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695]"; experiment "EXISTENCE:mutant phenotype:GO:0045944 positive regulation of transcription by RNA polymerase II [PMID:12586695]"; gbkey "CDS"; gene "PAU8"; locus_tag "YAL068C"; note "hypothetical protein; member of the seripauperin multigene family encoded mainly in subtelomeric regions"; product "seripauperin PAU8"; protein_id "NP_009332.1"; exon_number "1";
NC_001133.9 RefSeq  start_codon 2167    2169    .   -   0   gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "SGD:S000002142"; db_xref "GeneID:851229"; experiment "EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695]"; experiment "EXISTENCE:mutant phenotype:GO:0045944 positive regulation of transcription by RNA polymerase II [PMID:12586695]"; gbkey "CDS"; gene "PAU8"; locus_tag "YAL068C"; note "hypothetical protein; member of the seripauperin multigene family encoded mainly in subtelomeric regions"; product "seripauperin PAU8"; protein_id "NP_009332.1"; exon_number "1";
NC_001133.9 RefSeq  stop_codon  1807    1809    .   -   0   gene_id "YAL068C"; transcript_id "NM_001180043.1"; db_xref "SGD:S000002142"; db_xref "GeneID:851229"; experiment "EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695]"; experiment "EXISTENCE:mutant phenotype:GO:0045944 positive regulation of transcription by RNA polymerase II [PMID:12586695]"; gbkey "CDS"; gene "PAU8"; locus_tag "YAL068C"; note "hypothetical protein; member of the seripauperin multigene family encoded mainly in subtelomeric regions"; product "seripauperin PAU8"; protein_id "NP_009332.1"; exon_number "1";
```

[Table of contents](#table-of-contents)

---

### FASTQ format

The final input file format we will describe is the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format. The FASTQ format is similar to the FASTA format described in the previous section, but contains additional information describing the quality of the sequence:

```
cd /fastscratch/$USER/rnaseq/data

$ head SRR11849386.fastq
@SRR11849386.1 1 length=75
GGTCCNCCTCCCTAACGGGAGGGGGTCCNTCACTCCTTTCTTAAAATAAATTATTATTAATAATTATTATATAAT
+SRR11849386.1 1 length=75
AAAAA#EEEEEEAA6EEEEEEEEEEEAA#A/AEEEEAEEEEEEEEEEEEEEEEEEEEEEAEEAEE/AEEEEE6AE
@SRR11849386.2 2 length=74
CTCTTNATCAAATGGGTGGTGTAACCAGNAATCTTGTTTCTCAATCTCTTGGATTGGATAGTGGCGATTTCATC
+SRR11849386.2 2 length=74
AAAAA#AEEEEEEEEEEEEEEEEEEEAE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEE
@SRR11849386.3 3 length=75
TTATANATGTCCTCTATTCTAATACATACTTTTTTACTTTTGAAATAACAAAATATTCTTTTTATTTTTGCGAGT
```

In the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format, each sequence consists of exactly 4 lines. The first line of a record is the description line, which always has the character '`@`' in the first position. Just like the FASTA format, the first word (space delimited) on the description line corresponds to the sequence identifier. Subsequent words on this line are optional and fill out the sequence description. The second line is the sequence itself, which is never wrapped (always takes up exactly one line). The third (arguably pointless) line has a '`+`' character in the first position, and optionally reiterates the description line contents. The fourth line sets the FASTQ format apart from the FASTA format, by including sequence quality information. The sequence quality is encoded as a string of characters which is the same length as the sequence listed on the second line of the record. Each character on the fourth line encodes the [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score), or `Q`, which encodes an estimate of the probability `P`, that the corresponding base in the sequence has been identified (or 'called') incorrectly:

```
## using R syntax:
Q = -10 * log10(P)
```

or, conversely:

```
## more R syntax:
P = 10 ^ (-Q / 10)
```

The Phred score for a given sequence residue is encoded as a single printable character that appears in the corresponding position of the quality string. The character set used in this context is called the [ASCII](https://ss64.com/ascii.html) character set. This set includes both printable (e.g. `A`, `a`, `1`, `0`, `!`, `=`, etc.) and non-printable (e.g. space, backspace, tab, carriage return, etc.) characters. If you look at the table found by following the link above, you will see that there is a contiguous stretch consisting of only printable characters extending from `!` (corresponding to the decimal integer 33) to `~` which corresponds to the decimal integer 126). In order to map the Phred scores into this range, we add 33 to the rounded Phred score and print the corresponding character. We can see this in action by using a little R to decode the values. For the first sequence in the FASTA file shown above, the estimated quality codes for the first 7 bases `GGTCCNC` are `AAAAA#E`:

```
## a bit of R:

## convert ASCII character into decimal value:

> as.integer(charToRaw('A'))
[1] 65

## convert ASCII character into Phred score Q; the outer parentheses print value assigned:

> (Q = as.integer(charToRaw('A')) - 33)
[1] 32

## convert Phred score Q into probability of bad call P:

> 10 ^ (-Q / 10)
[1] 0.0006309573

## put it all together into a function that takes an ascii character and spits out P:

f <- function(ascii) {
  Q = as.integer(charToRaw(ascii)) - 33
  10 ^ (-Q / 10)
}

## our quality string started with AAAAA#E:

> f('A')
[1] 0.0006309573         ## about 0.06% chance bad call (approx. 1 in 1500)

> f('#')
[1] 0.6309573            ## >50% chance bad call -- no wonder called 'N'

> f('E')
[1] 0.0002511886         ## about 0.02% chance bad call (approx. 1 in 5000)

## the first and last characters in the allowable range:

> f('!')
[1] 1                    ## worst score possible; probability of bad call is 1 !!!

> f('~')
[1] 5.011872e-10         ## best score possible
```

[Table of contents](#table-of-contents)

---

Next: [Check sequence quality](fastqc.md)  
Previous: [Challenges to DE analysis](challenges.md)  

