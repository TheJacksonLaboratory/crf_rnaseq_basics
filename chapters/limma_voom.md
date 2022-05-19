# Lima-voom

## WORK IN PROGRESS

---

### Table of contents

- [Introduction](#introduction)
- [Algorithm details](#algorithm-details)
- [How to](#how-to)

#### Introduction

We will next do differential expression analysis using the R package limma and the 'voom' approach, whose distinctive innovation is described [here](https://doi.org/10.1186/gb-2014-15-2-r29).

---

#### Algorithm details

---

#### How to

```
rm(list=ls())

setwd("C:/Users/user15/data/rnaseq")

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

> table(i.keep, useNA='ifany')
FALSE  TRUE 
38518 22148 

dat <- dat[i.keep, , keep.lib.sizes=F]
dat <- edgeR::calcNormFactors(dat)
des <- model.matrix(~meta$disease) 

## mean vs. variance plot:
v <- limma::voom(dat, design=des, plot=T)

fit <- limma::lmFit(v, design=des)
fit <- limma::eBayes(fit)

## for 1 coef: moderated t-test; for >1 coef: moderated F-test:

> limma::topTable(fit, coef=1, number=5) 
                   logFC  AveExpr        t      P.Value    adj.P.Val        B
ENSG00000075624 14.24467 14.08730 61.30493 9.980136e-47 1.105200e-42 94.74364
ENSG00000204592 12.69477 12.57254 61.37031 9.495292e-47 1.105200e-42 94.69762
ENSG00000210082 17.27324 17.34214 60.32781 2.113970e-46 1.560673e-42 94.15766
ENSG00000166710 14.82770 14.82363 59.34913 4.537087e-46 2.512185e-42 93.42557
ENSG00000087086 11.69968 11.55681 58.98954 6.025411e-46 2.669016e-42 92.99029

> limma::topTable(fit, coef=2, number=5)
                    logFC   AveExpr         t     P.Value adj.P.Val         B
ENSG00000071553 -1.459662  5.480389 -3.241886 0.002178705 0.9997236 -4.461950
ENSG00000075624 -1.133079 14.087303 -2.281022 0.027093744 0.9997236 -4.472515
ENSG00000182551 -1.291163  6.773107 -2.585222 0.012873379 0.9997236 -4.489661
ENSG00000140396 -1.844712  5.155285 -3.005510 0.004233042 0.9997236 -4.495302
ENSG00000167996 -1.146207 12.223549 -2.115996 0.039636790 0.9997236 -4.496916

> limma::topTable(fit, coef=3, number=5)
                   logFC  AveExpr        t      P.Value adj.P.Val         B
ENSG00000127951 2.654252 3.922963 3.628703 0.0006970853 0.6453993 -4.426928
ENSG00000139318 2.039336 4.104504 3.402679 0.0013668992 0.6453993 -4.446585
ENSG00000179344 1.794213 3.563195 3.628161 0.0006982295 0.6453993 -4.446663
ENSG00000163131 1.905036 6.512962 2.898385 0.0056695075 0.6453993 -4.453638
ENSG00000182718 1.507547 5.057839 3.117518 0.0031001287 0.6453993 -4.455375

> limma::topTable(fit, coef=2:3, number=5)
                meta.diseasemucinousadenocarcinoma meta.diseaseserousadenocarcinoma     AveExpr         F      P.Value adj.P.Val
ENSG00000053918                         -2.1854512                         2.493383  0.61068439 12.921674 3.332356e-05 0.7380502
ENSG00000120738                         -3.1621053                         1.910497  0.86931824  9.279837 3.982218e-04 0.7707887
ENSG00000214063                         -3.1963715                         0.603909  0.71337640  9.106378 4.512020e-04 0.7707887
ENSG00000197168                          2.9678828                        -1.766876 -0.03622534  9.051885 4.693229e-04 0.7707887
ENSG00000128268                         -0.8406509                        -3.132850 -0.87368405  8.904265 5.223256e-04 0.7707887

```

[Table of contents](#table-of-contents)

---
