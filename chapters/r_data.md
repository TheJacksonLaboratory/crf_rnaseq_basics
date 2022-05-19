## Import data into R

- [Home](../README.md)
- [Experimental Design](design.md)
- [Challenges to DE analysis](challenges.md)
- [Inputs](inputs.md)
- [Check sequence quality](fastqc.md)
- [Map reads](mapping.md)
- [Generate expression matrix](count_matrix.md)
- **[Import data into R](r_data.md)** *(You are here)*
- [Screen for outlier samples](outliers.md)
- [edgeR](edger.md)
- [DESeq2](deseq2.md)

---

Once we have our feature count matrix (from the previous chapter), we can proceed with the rest of the analysis on the JAX HPC cluster or on a laptop. Put `feature_counts.txt` (from `/fastscratch/$USER/rnaseq/feature_counts/`) and `metadata.GSE151185.txt` (from `/fastscratch/$USER/rnaseq/data/`) into a new working folder (perhaps on your local computer) and launch R. 

On the JAX HPC cluster you can do the following:

```
## set up the working directory and link in the input data:

cd /fastscratch/$USER/rnaseq
mkdir R
ln -s /fastscratch/$USER/rnaseq/feature_counts/feature_counts.txt R/
ln -s /fastscratch/$USER/rnaseq/data/metadata.GSE151185.txt R/
cd R

## how to get R going (with all the packages you will need) on the JAX HPC cluster:

module load singularity
singularity exec /projects/researchit/crf/containers/crf_rnaseq.sif R
```

On your laptop, download `feature_counts.txt` and `metadata.GSE151185.txt` into an empty working directory, change into that directory and launch R. Once in R, you may need to change to that working directory using the `setwd("c:/path/to/my/data_directory")` R command. Remember to use the forward slash `/` as your path delimiter in R, even on Windows. You can get help for any R command `cmd()` by entering `?cmd`, for instance, to see how `setwd()` works, you can enter `?setwd`.

Once we have R going on either our laptop or on the JAX HPC cluster, we can import our data and begin to explore it:  

```
## R: comments start with '#':

## Read in the metadata file; include header row; look for tab separator, 
##   do not convert character strings to factors (as.is=T):

meta <- read.table("metadata.GSE151185.txt", header=T, sep="\t", as.is=T)

## Read in the counts matrix:

dat <- read.table("feature_counts.txt", header=T, sep="\t", as.is=T)

## how much RAM do our new objects consume:

> object.size(meta)
9680 bytes

> object.size(dat)
2084384 bytes
```

Let us take a look at `meta`, to make sure everything looks ok and to learn more about the variables we have to work with:

```
## the meta object is a 'data.frame', which is a rectangular table with
##   rows corresponding to observations (samples) and columns corresponding 
##   to variables:

> class(meta)
[1] "data.frame"

## this table has 27 rows and 7 columns:

> dim(meta)
[1] 27  7

## let's look at the top 6 (default for head()) rows:

> head(meta)
      Sample        SRX         SRR      Label treatment time rep
1 GSM4568115 SRX8399630 SRR11849386 Waterlog_1     water    0   1
2 GSM4568116 SRX8399631 SRR11849387 Waterlog_2     water    0   2
3 GSM4568117 SRX8399632 SRR11849388 Waterlog_2     water    0   3
4 GSM4568118 SRX8399606 SRR11849389  Water24_1     water   24   1
5 GSM4568119 SRX8399607 SRR11849390  Water24_2     water   24   2
6 GSM4568120 SRX8399608 SRR11849391  Water24_3     water   24   3

## our variables of interest are 'treatment' and 'time':

## 9 replicates for each of three treatments (CR, NR, and water):

> table(meta$treatment)
   CR    NR water 
    9     9     9 

## 9 replicates for each of three timepoints (0, 24 and 96 hours):

> table(meta$time)
 0 24 96 
 9  9  9 

## the experiment is 'balanced': there are the same number of  
##   replicates (n=3) for each variable combination, which is the 
##   ideal experimental layout (though replication is low!!!); we 
##   want at least n=5, and more ideally n>=10:

> table(meta$treatment, meta$time)
        0 24 96
  CR    3  3  3
  NR    3  3  3
  water 3  3  3
```

Next we will look at the count matrix. The rows of this matrix correspond to genes, and the columns correspond to samples. This matrix includes the `Geneid` as the first column, followed by columns specifying exon coordinates for the gene. These are followed by one column for each sample, where each cell contains the number reads mapped to the corresponding row's gene in this column's sample:

```
> class(dat)
[1] "data.frame"

> dim(dat)
[1] 6458   33

## how the columns are labeled; we just need Geneid and the counts for each sample,
##   not the columns with the gene's coordinates (columns 2-6):

> names(dat)
 [1] "Geneid"                     "Chr"                       
 [3] "Start"                      "End"                       
 [5] "Strand"                     "Length"                    
 [7] "SRR11849386Aligned.out.bam" "SRR11849387Aligned.out.bam"
 [9] "SRR11849388Aligned.out.bam" "SRR11849389Aligned.out.bam"
[11] "SRR11849390Aligned.out.bam" "SRR11849391Aligned.out.bam"
[13] "SRR11849392Aligned.out.bam" "SRR11849393Aligned.out.bam"
[15] "SRR11849394Aligned.out.bam" "SRR11849395Aligned.out.bam"
[17] "SRR11849396Aligned.out.bam" "SRR11849397Aligned.out.bam"
[19] "SRR11849398Aligned.out.bam" "SRR11849399Aligned.out.bam"
[21] "SRR11849400Aligned.out.bam" "SRR11849401Aligned.out.bam"
[23] "SRR11849402Aligned.out.bam" "SRR11849403Aligned.out.bam"
[25] "SRR11849404Aligned.out.bam" "SRR11849405Aligned.out.bam"
[27] "SRR11849406Aligned.out.bam" "SRR11849407Aligned.out.bam"
[29] "SRR11849408Aligned.out.bam" "SRR11849409Aligned.out.bam"
[31] "SRR11849410Aligned.out.bam" "SRR11849411Aligned.out.bam"
[33] "SRR11849412Aligned.out.bam"

## first five genes in first two samples (columns 7 and 8):

> dat[1:5, 1:8]
     Geneid         Chr Start   End Strand Length SRR11849386Aligned.out.bam SRR11849387Aligned.out.bam
1   YAL068C NC_001133.9  1807  2169      -    363                          5                          6
2 YAL067W-A NC_001133.9  2480  2707      +    228                          5                          5
3   YAL067C NC_001133.9  7235  9016      -   1782                         93                         75
4   YAL065C NC_001133.9 11565 11951      -    387                         24                         15
5 YAL064W-B NC_001133.9 12046 12426      +    381                         31                         27

## We need Geneid, but let us make it a 'rowname' rather than a column;
##   once assigned as a rowname, we will no longer need the Geneid column:

rownames(dat) <- dat$Geneid             ## set rownames of dat to Geneid column of dat

## now we can get rid of columns 1-6;
##   the ':' operator flanked by integers defines a range:

> 1:6
[1] 1 2 3 4 5 6

## the labels of columns we no longer need; data.frame column labels accessible with 
##   names() or colnames(); can only use colnames() for matrices:

> names(dat)[1:6]
[1] "Geneid" "Chr"    "Start"  "End"    "Strand" "Length"

## the colums we do need (use negative integers to omit columns):

> names(dat)[-(1:6)]
 [1] "SRR11849386Aligned.out.bam" "SRR11849387Aligned.out.bam"
 [3] "SRR11849388Aligned.out.bam" "SRR11849389Aligned.out.bam"
 [5] "SRR11849390Aligned.out.bam" "SRR11849391Aligned.out.bam"
 [7] "SRR11849392Aligned.out.bam" "SRR11849393Aligned.out.bam"
 [9] "SRR11849394Aligned.out.bam" "SRR11849395Aligned.out.bam"
[11] "SRR11849396Aligned.out.bam" "SRR11849397Aligned.out.bam"
[13] "SRR11849398Aligned.out.bam" "SRR11849399Aligned.out.bam"
[15] "SRR11849400Aligned.out.bam" "SRR11849401Aligned.out.bam"
[17] "SRR11849402Aligned.out.bam" "SRR11849403Aligned.out.bam"
[19] "SRR11849404Aligned.out.bam" "SRR11849405Aligned.out.bam"
[21] "SRR11849406Aligned.out.bam" "SRR11849407Aligned.out.bam"
[23] "SRR11849408Aligned.out.bam" "SRR11849409Aligned.out.bam"
[25] "SRR11849410Aligned.out.bam" "SRR11849411Aligned.out.bam"
[27] "SRR11849412Aligned.out.bam"

## get rid of the columns we do not need; remember indexing is dat[row, column];
##   dat[, -(1:6)] keeps all rows and drops columns 1-6:

dat <- dat[, -(1:6)]

> dim(dat)
[1] 6458   27

## now we can look at our expression matrix; too many columns for head();
##   will index instead to view rows 1-6 and columns 1-3; note that we can see the rownames
##   (Geneids) and colnames (BAM file names for each sample):

> dat[1:6, 1:3]
          SRR11849386Aligned.out.bam SRR11849387Aligned.out.bam SRR11849388Aligned.out.bam
YAL068C                            5                          6                          1
YAL067W-A                          5                          5                          3
YAL067C                           93                         75                        125
YAL065C                           24                         15                         28
YAL064W-B                         31                         27                         40
YAL064C-A                         87                         55                        109
```

We will need to be able to associate the rows (samples) of `meta` with the sample columns of `dat`. To make that easier, we will reformat the column labels of `dat` to match the `SRR` column of `meta`, by getting rid of the trailing `Aligned.out.bam` substring:

```
## current column 1-4 labels:

> names(dat)[1:4]
[1] "SRR11849386Aligned.out.bam" "SRR11849387Aligned.out.bam" "SRR11849388Aligned.out.bam" "SRR11849389Aligned.out.bam"

## substitute 'Aligned.out.bam' with ''; '$' at end of pattern ensures only matches at end of string:

> sub('Aligned.out.bam$', '', names(dat)[1:4])
[1] "SRR11849386" "SRR11849387" "SRR11849388" "SRR11849389"

## process all the column labels:

names(dat) <- sub('Aligned.out.bam$', '', names(dat))

## all look good:

> names(dat)
 [1] "SRR11849386" "SRR11849387" "SRR11849388" "SRR11849389" "SRR11849390"
 [6] "SRR11849391" "SRR11849392" "SRR11849393" "SRR11849394" "SRR11849395"
[11] "SRR11849396" "SRR11849397" "SRR11849398" "SRR11849399" "SRR11849400"
[16] "SRR11849401" "SRR11849402" "SRR11849403" "SRR11849404" "SRR11849405"
[21] "SRR11849406" "SRR11849407" "SRR11849408" "SRR11849409" "SRR11849410"
[26] "SRR11849411" "SRR11849412"

## what the table looks like. rownames on left; colnames on top; 
##   cells have integers:

> dat[1:8, 1:5]
          SRR11849386 SRR11849387 SRR11849388 SRR11849389 SRR11849390
YAL068C             5           6           1           8           4
YAL067W-A           5           5           3           6           1
YAL067C            93          75         125         477         437
YAL065C            24          15          28          28          31
YAL064W-B          31          27          40          36          36
YAL064C-A          87          55         109          98          84
YAL064W             7           5           3           8           6
YAL063C-A          18          30          15          64          48
```

We have prepared our data for analysis by getting rid of the unneeded columns from the count matrix and ensuring that the sample labels in that matrix can be easily cross-referenced to rows in the metadata file. This will allow us to integrate information from these two tables. In order to not have to repeat this process, we will save this intermediate result in the relatively human readable tab-delimited format (takes a bit more space, but more versatile):

```
write.table(dat, file="count_matrix.txt", sep="\t")
```

---

Next: [Screen for outlier samples](outliers.md)    
Previous: [Generate expression matrix](count_matrix.md)  
