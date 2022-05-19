# Lima-trend

## WORK IN PROGRESS

---

### Table of contents

- [Introduction](#introduction)
- [Algorithm details](#algorithm-details)
- [How to](#how-to)

## Introduction

We will next do differential expression analysis with the R package [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) using the 'trend' approach. The limma package is older than the other packages described here. Limma was one of the most commonly used tools for differential expression analysis of data generated using [cDNA microarrays](https://en.wikipedia.org/wiki/DNA_microarray), and was later adapted to use for RNA-seq. Limma provided one of the first effective approaches to moderating effects of dispersion estimate imprecision when testing many genes using small sample sizes (discussed in the [Challenges](#challenges-to-de-analysis) section), that were later incorporated into the edgeR package. The creator of limma (Gordon Smyth) was a coathor of the [edgeR](#edger) paper and software, so it may come as no surprise that the two packages are often used together, with edgeR used for data import and inter-sample normalization, followed by use of limma for model fitting and hypothesis testing. Limma appears less popular than either edgeR or DESeq2, but is very fast and can prove useful in cases where cross-validation, empirical bootstrapping, or jackkniffing are conducted.

---

## Algorithm details

When analysing microarrays, it is quite common to use [quantile normalization](https://doi.org/10.1093/bioinformatics/19.2.185) for inter-sample normalization. This method is also sometimes used for normalization of RNA-seq data. However, the method does not work as well for RNA-seq data as other methods discussed in this tutorial. The reason is likely related to the fact that different genes compete for a fixed size read pool in RNA-seq experiments, so if one gene becomes more highly expressed, it uses up more of the available reads, leaving fewer reads available for other genes, causing their read counts to decline. This introduces a technical (rather than biological) dependency between different gene signals that is not present in microarray experiments. In microarray experiments, each gene/isoform has its own spot (or set of spots), and binding of flourescently labeled RNA to one spot does not appreciably affect measurements at other spots. Therefore, we do not recommend use of quantile normalization for RNA-seq data, and in the example below, we use the edgeR package to import the data and perform inter-sample normalization, then pass the resulting object directly to limma.

 The distinctive part of the limma package is its use of weighted linear regression for model fitting. 

---

## How to

```
setwd("C:/Users/$USER/data/rnaseq")

x <- read.table("feature_counts.txt", skip=1, header=T, sep="\t", as.is=T)
genes <- x[, 1:6]
rownames(x) <- x$Geneid
x <- x[, 7:ncol(x)]
names(x) <- sub("Aligned.out.bam", "", names(x))

meta <- read.table("metadata.txt", header=T, sep="\t", as.is=T)
rownames(meta) <- meta$sra_id
meta <- meta[names(x), ]

i <- meta$disease %in% c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma')
meta <- meta[i, ]
rownames(meta) <- NULL

x <- x[, meta$sra_id]

## ensure that control condition is first factor level:

meta$disease <- factor(meta$disease, 
  levels=c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma'))

dat <- edgeR::DGEList(counts=x)
min.samples <- 3
min.counts <- 10
i.keep <- apply(x, 1, function(v) sum(v >= min.counts, na.rm=T) >= min.samples)

dat <- dat[i.keep, , keep.lib.sizes=F]
dat <- edgeR::calcNormFactors(dat)

log.cpm <- edgeR::cpm(dat, log=T, prior.count=3)
des <- model.matrix(~meta$disease) 

fit <- limma::lmFit(log.cpm, design=des)
fit <- limma::eBayes(fit, trend=TRUE)

## for 1 coef: moderated t-test; for >1 coef: moderated F-test:

> limma::topTable(fit, coef=1, number=5) 
                   logFC  AveExpr        t      P.Value    adj.P.Val         B
ENSG00000166710 14.82690 14.82366 151.2523 7.881485e-63 1.745591e-58 103.52432
ENSG00000210082 17.27087 17.34214 144.4361 6.319143e-62 6.997818e-58 103.00644
ENSG00000211459 14.68722 14.77767 120.4374 2.297503e-58 1.696170e-54 100.59242
ENSG00000198888 15.38447 15.46347 113.0165 4.041692e-57 2.237885e-53  99.59527
ENSG00000198727 13.36827 13.46238 110.6800 1.036383e-56 4.084926e-53  99.24942

> limma::topTable(fit, coef=2, number=5)
                    logFC   AveExpr         t     P.Value adj.P.Val         B
ENSG00000176200  2.196273 0.3606955  3.400193 0.001416266 0.9998574 -4.506011
ENSG00000005882 -2.268117 2.2841232 -3.290980 0.001940582 0.9998574 -4.511111
ENSG00000145982 -1.483557 3.5725335 -3.245564 0.002208909 0.9998574 -4.513221
ENSG00000256618  2.452429 1.3115048  3.216633 0.002397766 0.9998574 -4.514561
ENSG00000137878 -1.363261 6.0449029 -3.168569 0.002745678 0.9998574 -4.516780

> limma::topTable(fit, coef=3, number=5)
                    logFC    AveExpr         t     P.Value adj.P.Val         B
ENSG00000204710 -2.524127 3.22903165 -3.394811 0.001438586 0.9005333 -4.518016
ENSG00000187714  1.483744 0.03530998  3.011060 0.004249252 0.9005333 -4.533387
ENSG00000256618  2.450653 1.31150478  2.990140 0.004499121 0.9005333 -4.534209
ENSG00000279267 -2.150084 2.07019454 -2.700728 0.009701401 0.9005333 -4.545356
ENSG00000166311 -1.391009 3.79976215 -2.687312 0.010042778 0.9005333 -4.545861

> limma::topTable(fit, coef=2:3, number=5)
                meta.diseasemucinousadenocarcinoma meta.diseaseserousadenocarcinoma   AveExpr        F     P.Value adj.P.Val
ENSG00000256618                         2.45242875                        2.4506527 1.3115048 7.997283 0.001061944 0.9710095
ENSG00000176200                         2.19627325                       -0.2654143 0.3606955 6.395720 0.003580032 0.9710095
ENSG00000204710                         0.04046983                       -2.5241275 3.2290316 6.064302 0.004642337 0.9710095
ENSG00000005882                        -2.26811693                        0.2100621 2.2841232 5.900531 0.005284262 0.9710095
ENSG00000095637                         1.88740074                       -1.4410516 2.3588357 5.878838 0.005375995 0.9710095
```

[Table of contents](#table-of-contents)

---
