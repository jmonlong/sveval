# sveval

[![Build Status](https://travis-ci.com/jmonlong/sveval.svg?branch=master)](https://travis-ci.com/jmonlong/sveval)

Functions to compare a SV call sets against a truth set.

## Installation

```r
install.packages('devtools') ## if devtools not already installed
source('http://bioconductor.org/biocLite.R')
biocLite('jmonlong/sveval')
```

Or to install locally (e.g. in a HPC)

```r
.libPaths('~/R/library/')
install.packages('devtools') ## if devtools not already installed
source('http://bioconductor.org/biocLite.R')
biocLite('jmonlong/sveval')
```

## Usage

```r
library(sveval)
eval.o = svevalOl('calls.vcf', 'truth.vcf')
eval.o$eval # data.frame with results using all variants
plot_prcurve(eval.o$curve)
```

Outputs a list with a data.frame with TP, FP, TN, precision, recall and F1 for all variants and for each SV type, and a another data.frame with the results using increasing quality thresholds to make a precision-recall curve.

Some of the most important other parameters:

- `max.ins.dist=` maximum distance for insertions to be clustered. Default is 20.
- `min.cov=` the minimum coverage to be considered a match. Default is 0.5
- `min.del.rol=` minimum reciprocal overlap for deletions. Default is 0.1
- `min.size=` the minimum SV size to be considered. Default 0.
- `bed.regions=` If non-NULL, a GRanges object or path to a BED file (no headers) with regions of interest.
- `outfile=` the TSV file to output the results. If NULL (default), returns a data.frame.
- `ins.seq.comp=TRUE` compare sequence instead of insertion sizes. Default is *FALSE*.
- `check.inv` should the sequence of MNV be compared to identify inversions. Default is *FALSE*.
- `geno.eval`/`merge.hets`/`stitch.hets` options for genotype evaluation, see below.

See full list of parameters in the [manual](docs/sveval-manual.pdf) or by typing `?svevalOl` in R.

### Precision-recall curve comparing multiple methods

```r
eval.1 = svevalOl('calls1.vcf', 'truth.vcf')
eval.2 = svevalOl('calls2.vcf', 'truth.vcf')
plot_prcurve(list(eval.1$curve, eval.2$curve), labels=c('method1', 'method2'))
```

Or if the results were written in files:

```r
plot_prcurve(c('methods1-prcurve.tsv', 'methods2-prcurve.tsv'), labels=c('method1', 'method2'))
```

### Genotype evaluation

By default sveval doesn't take the genotype into account, more a calling evaluation than a genotype evaluation.
To compare genotype, the evaluation can be performed separately for heterozygous and homozygous variants.
Before doing that it sometimes help to merge very similar hets into homs.
To a lower extent, it also helps to stitch fragmented hets before trying to merge them into homs.
The relevant parameters in `svevalOl` are:

- `geno.eval=TRUE` compare hets/homs separately.
- `stitch.hets=TRUE` stitch fragmented hets.
- `stitch.dist` the maximum distance between two hets to be stitched. Default 20 bp.
- `merge.hets=TRUE` merge hets into hom before comparison.
- `merge.rol` the minimum reciprocal overlap between two hets to be merged. Default is 0.8.

Hence, the **recommended command for genotype evaluation**:

```r
eval.o = svevalOl('calls.vcf', 'truth.vcf', geno.eval=TRUE, stitch.hets=TRUE, merge.hets=TRUE)
```

## Methods

- For deletions, at least 50% coverage and at least 10% reciprocal overlap.
- For insertions, size of nearby insertions (+- 20 bp) at least as much as 50% the size of insertion. Or comparing inserted sequence (sequence similarity instead of size).
- For inversions, same as deletions. If using REF/ALT sequences (i.e. not symbolic ALT), inversions are variants longer than 10 bp where the reverse complement of ALT matches REF at least 80%.

![](docs/ol-cartoon.svg)

## Docker

A docker image of R with this package installed is available [here](https://hub.docker.com/r/jmonlong/sveval/).
