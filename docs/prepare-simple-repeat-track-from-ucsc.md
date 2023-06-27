# Prepare simple repeat track

``` r
library(GenomicRanges)
library(dplyr)
```

## Download simple repeat track from the UCSC Genome Browser FTP

This annotation originally came from TRF.

``` r
GENOME = 'GRCh38'
## url of the file on the UCSC Genome Browser repository
sr.url = NULL
if(GENOME=='GRCh38'){
  sr.url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz'
} else if (GENOME=='GRCh37') {
  sr.url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz'
}

if(!file.exists(paste0('simpleRepeat_', GENOME, '.txt.gz'))){
  download.file(sr.url, paste0('simpleRepeat_', GENOME, '.txt.gz'))
}
```

## Read and keep coordinate columns only

``` r
sr = read.table(paste0("simpleRepeat_", GENOME, ".txt.gz"), as.is = TRUE, sep = "\t")
sr = with(sr, GRanges(V2, IRanges(V3, V4)))
sr
```

    ## GRanges object with 1049715 ranges and 0 metadata columns:
    ##                        seqnames        ranges strand
    ##                           <Rle>     <IRanges>  <Rle>
    ##         [1]                chr1   10000-10468      *
    ##         [2]                chr1   10627-10800      *
    ##         [3]                chr1   10757-10997      *
    ##         [4]                chr1   11225-11447      *
    ##         [5]                chr1   11271-11448      *
    ##         ...                 ...           ...    ...
    ##   [1049711] chrY_MU273398v1_fix 846718-846801      *
    ##   [1049712] chrY_MU273398v1_fix 852183-852214      *
    ##   [1049713] chrY_MU273398v1_fix 855179-855214      *
    ##   [1049714] chrY_MU273398v1_fix 855212-855261      *
    ##   [1049715] chrY_MU273398v1_fix 861797-861856      *
    ##   -------
    ##   seqinfo: 702 sequences from an unspecified genome; no seqlengths

Total of annotated regions: 345.99Mbp. However, we can see that some
annotated regions overlap significantly.

## Reduce overlapping repeats

As shown above in the first annotated regions, the same region might be
annotated multiple times. We “merge” them to get larger annotated
regions. Larger regions mean that sveval will be able to move SVs more
in those regions when attempting to match them.

``` r
sr = reduce(sr)
sr
```

    ## GRanges object with 702648 ranges and 0 metadata columns:
    ##                       seqnames        ranges strand
    ##                          <Rle>     <IRanges>  <Rle>
    ##        [1]                chr1   10000-10468      *
    ##        [2]                chr1   10627-10997      *
    ##        [3]                chr1   11225-11448      *
    ##        [4]                chr1   19305-19443      *
    ##        [5]                chr1   20828-20863      *
    ##        ...                 ...           ...    ...
    ##   [702644] chrY_MU273398v1_fix 846429-846466      *
    ##   [702645] chrY_MU273398v1_fix 846718-846801      *
    ##   [702646] chrY_MU273398v1_fix 852183-852214      *
    ##   [702647] chrY_MU273398v1_fix 855179-855261      *
    ##   [702648] chrY_MU273398v1_fix 861797-861856      *
    ##   -------
    ##   seqinfo: 702 sequences from an unspecified genome; no seqlengths

Total of annotated regions: 152.1Mbp.

## Write gzipped bed

Keep only primary sequences (no `_` in the sequence name).

``` r
outbed = gzfile(paste0("simpleRepeat_", GENOME, ".bed.gz"), "w")
sr %>%
    as.data.frame %>%
    select(seqnames, start, end) %>%
    filter(!grepl("_", seqnames)) %>%
    mutate(seqnames = ifelse(GENOME == "hg19", gsub("chr", "", seqnames), seqnames)) %>%
    write.table(file = outbed, col.names = FALSE, quote = FALSE, row.names = FALSE,
        sep = "\t")
close(outbed)
```
