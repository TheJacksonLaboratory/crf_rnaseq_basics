## Generate expression matrix

- [Home](../README.md)
- [Experimental Design](design.md)
- [Challenges to DE analysis](challenges.md)
- [Inputs](inputs.md)
- [Check sequence quality](fastqc.md)
- [Map reads](mapping.md)
- **[Generate expression matrix](count_matrix.md)** *(You are here)*
- [Import data into R](r_data.md)
- [Screen for outlier samples](outliers.md)
- [edgeR](edger.md)
- [DESeq2](deseq2.md)

---

The next step in our analysis will be the generation of a gene expression matrix, with each gene represented by a row and each sample represented by a column. Each cell should then contain an integer (greater than or equal to zero) corresponding to the number of reads from the column's sample that were mapped to the row's gene. Additional columns can be included with this basic metadata, containing information such as the identifier for each gene. The resulting matrix file will be orders of magnitude smaller than the BAM formatted alignment files on which it is based. Once the gene expression matrix file is generated, it can usually be transferred to a laptop, with the rest of the analysis taking place there.

We will generate an expression matrix using the [featureCounts](http://subread.sourceforge.net/featureCounts.html) program from the [Subreads](http://subread.sourceforge.net/) R package; however `featureCounts` is available as a stand-alone program (in our container), so we do not need to invoke R to run it. The basic command is:

```
featureCounts \
  -T $cpus \
  -Q $mapq_min \
  --minOverlap $min_overlap \
  --primary \
  --ignoreDup \
  -a "$gtf" \
  -o "$out_file" \
  $bam1 $bam2 $bam3 ...
```

Here `$cpus` specifies the number of compute threads to use, and `$mapq_min` is the minimum acceptable MAPQ score required of an alignment for it to be counted (we recommend a cutoff of `10`, per the discussion at the end of the [Map reads](mapping.md) chapter). By the default settings, if a read overlaps an exon by at least 1 residue, the read is counted for that exon's gene. This stringency can be controlled by the `--minOverlap` parameter, which we typically set to `10`. The `--primary` flag filters out alignments whose flag (the second column of the SAM format) indicates a secondary alignment (the 0x100 bit is set). The `--ignoreDup` flag causes the program to ignore reads flagged as duplicates (the 0x400 bit is set), which is intended to filter out reads from sequences that were artifactually duplicated due to PCR steps during library preparation or due to imaging artifacts sometimes encountered during sequencing. You may not want to use the `--ignoreDup` flag if your sequencing is deep enough relative to by chance create many exactly matching reads of highly expressed transcripts. If you do want to use it, you may have to use something like `samtools markdup` to mark (or just remove) reads deemed 'duplicates'. Samtools judges two reads duplicates if the 5'-end mapping coordinate of two reads coincides. In that case, the one with higher mapping quality is kept, while the other is marked as a duplicate or removed, depending on other flags you can select. The command is documented [here](http://www.htslib.org/doc/samtools-markdup.html) and is available in our container. The reference genome exons are provided as the GTF2.2 formatted file specified by `$gtf`, and the output file will be named `$out_file`. Finally, we provide our alignment BAM files as a sequence of filenames separated by spaces.

We'll parallelize the process using the convenience script `feature_counts.sh`. As usual, we'll first create a working directory for this step and link our inputs (the BAM files STAR alignment resulted in) into this directory:

```
mkdir -p /fastscratch/$USER/rnaseq/feature_counts
cd /fastscratch/$USER/rnaseq
ln -s $PWD/star/*.bam feature_counts
ln -s $PWD/reference/GCF_000146045.2_R64_genomic.gtf feature_counts
cp ~/opt/crf_rnaseq_basics/scripts/feature_counts.sh feature_counts
cd feature_counts

$ ls -l
total 0
-rw-r--r-- 1 user15 group12 1090 Feb 11 13:24 feature_counts.sh
lrwxrwxrwx 1 user15 group12   68 Feb 11 13:24 GCF_000146045.2_R64_genomic.gtf -> /fastscratch/user15/rnaseq/reference/GCF_000146045.2_R64_genomic.gtf
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849386Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849386Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849387Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849387Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849388Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849388Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849389Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849389Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849390Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849390Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849391Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849391Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849392Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849392Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849393Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849393Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849394Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849394Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849395Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849395Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849396Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849396Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849397Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849397Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849398Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849398Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849399Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849399Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849400Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849400Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849401Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849401Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849402Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849402Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849403Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849403Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849404Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849404Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849405Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849405Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849406Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849406Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849407Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849407Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849408Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849408Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849409Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849409Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849410Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849410Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849411Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849411Aligned.out.bam
lrwxrwxrwx 1 user15 group12   58 Feb 11 13:24 SRR11849412Aligned.out.bam -> /fastscratch/user15/rnaseq/star/SRR11849412Aligned.out.bam
```

The script `feature_counts.sh` looks for BAM files (any filename ending in '.bam') in the directory within which the command is invoked:

```
cd /fastscratch/$USER/rnaseq/feature_counts
sbatch feature_counts.sh GCF_000146045.2_R64_genomic.gtf feature_counts.txt
```

In addition to examining the usual `feature_counts.*.out` and `feature_counts.*.err` files to ensure the run completed successfully, you can examine the text file `feature_counts.txt.summary` (using `cat` or `less`) to see a summary of how many reads were successfully assigned for each sample:

```
$ cat feature_counts.txt.summary
Status  SRR11849386Aligned.out.bam  SRR11849387Aligned.out.bam  SRR11849388Aligned.out.bam  SRR11849389Aligned.out.bam  SRR11849390Aligned.out.bam   SRR11849391Aligned.out.bam  SRR11849392Aligned.out.bam  SRR11849393Aligned.out.bam  SRR11849394Aligned.out.bam   SRR11849395Aligned.out.bam  SRR11849396Aligned.out.bam  SRR11849397Aligned.out.bam  SRR11849398Aligned.out.bam  SRR11849399Aligned.out.bam   SRR11849400Aligned.out.bam  SRR11849401Aligned.out.bam  SRR11849402Aligned.out.bam  SRR11849403Aligned.out.bam   SRR11849404Aligned.out.bam  SRR11849405Aligned.out.bam  SRR11849406Aligned.out.bam  SRR11849407Aligned.out.bam   SRR11849408Aligned.out.bam  SRR11849409Aligned.out.bam  SRR11849410Aligned.out.bam  SRR11849411Aligned.out.bam  SRR11849412Aligned.out.bam
Assigned    16342681    17497269    35838437    13845100    11369505    18161039    20036174    10071770    18115361    12621037 10483651    6712267 9288604 7126399 104104334   7070585 31322046    11548162    7684067 13487088    6376777 10264440    11135794 8401709 7131917 6789650 9429856
Unassigned_Unmapped 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
Unassigned_Read_Type    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   00
Unassigned_Singleton    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   00
Unassigned_MappingQuality   3196440 3106546 9296951 2634695 2056317 3916434 3073717 1757000 2892032 2403894 1917942 6283939 2018903  1194559 20413953    1781363 5540363 2042922 1472641 2624875 3469153 2150824 2092632 1920063 1599316 1912963 1506001
Unassigned_Chimera  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
Unassigned_FragmentLength   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   00   0
Unassigned_Duplicate    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   00
Unassigned_MultiMapping 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   00
Unassigned_Secondary    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   00
Unassigned_NonSplit 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
Unassigned_NoFeatures   1166218 1504171 2147972 1503626 1362781 2196018 2462830 1607507 2434040 796435  751826  397052  992944   817113  10931841    710881  3191445 1199828 500727  794100  407871  902116  1083559 823222  616270  701102  902431
Unassigned_Overlapping_Length   57882   74704   116883  61735   54992   88220   93619   58514   95554   43066   40448   2119739996   32749   436623  29513   138830  49715   27099   42980   21470   37559   44434   34109   28824   29609   38390
Unassigned_Ambiguity    47235   52394   114297  66043   55963   77897   83200   47998   77148   37488   30772   25528   4558935888   532276  57052   195356  70737   22827   38246   23168   48456   51709   39099   29574   30004   38894
```

Creation of the feature count matrix summarizes many gigabytes worth of binary-formatted mapping data into a few dozen megabyte human-readable text file. You will only need this file and the metadata file for the experiment for subsequent steps, and can download it to your laptop. We will take a look at the format of the feature count matrix file in the [next chapter](r_data.md).

---

Next: [Import data into R](r_data.md)  
Previous: [Map reads](mapping.md)  
