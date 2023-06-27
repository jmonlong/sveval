## Wrapper around the svevalOl function in the sveval package.
## Arguments
## 1: calls vcf
## 2: truth vcf
## 3: sample name
## 4: look for inversions?
## 5: regions of interest (NA means whole genome)
## 6: simple repeat track (NA means disabled)
## 7: genotype evaluation?
## 8: output Rdata file
## 9: minimum overlap to match variants (optional, default:0.5)
args = commandArgs(TRUE)

## If specified, the minimum coverage/reciprocal overlap to match variants
if(length(args)==8){
  args = c(args, '0.5')
}

library(sveval)

## Regions of interest or whole-genome?
bed = NULL
if(args[5] != 'NA'){
  bed = args[5]
}

## simple repeat track
sr = NULL
if(args[6] != 'NA'){
  library(GenomicRanges)
  sr = read.table(args[6], as.is=TRUE, sep='\t')
  sr = GRanges(sr[,1], IRanges(sr[,2], sr[,3]))
}

## evaluation
eval.o = svevalOl(args[1], args[2], sample.name=args[3], check.inv=as.logical(args[4]),
                  max.ins.dist=100,
                  bed.regions=bed, min.ol=as.numeric(args[9]), 
                  geno.eval=as.logical(args[7]=='geno'),
                  stitch.hets=as.logical(args[7]=='geno'),
                  merge.hets=as.logical(args[7]=='geno'),
                  simprep=sr, min.size=50)

## save RData object
save(eval.o, file=args[8])

