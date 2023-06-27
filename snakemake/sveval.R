## Wrapper around the svevalOl function in the sveval package.
## Arguments
## 1: calls vcf
## 2: truth vcf
## 3: sample name
## 4: look for inversions?
## 5: regions of interest (NA means whole genome)
## 6: genotype evaluation?
## 7: output Rdata file
## 8: minimum overlap to match variants (optional, default:0.5)
args = commandArgs(TRUE)

## If specified, the minimum coverage/reciprocal overlap to match variants
if(length(args)==7){
  args = c(args, '0.5')
}

library(sveval)

## Regions of interest or whole-genome?
bed = NULL
if(args[5] != 'NA'){
  bed = args[5]
}

## evaluation
eval.o = svevalOl(args[1], args[2], sample.name=args[3], check.inv=as.logical(args[4]),
                  max.ins.dist=100,
                  bed.regions=bed, min.ol=as.numeric(args[8]), 
                  geno.eval=as.logical(args[6]=='geno'), stitch.hets=as.logical(args[6]=='geno'),
                  merge.hets=as.logical(args[6]=='geno'), min.size=50)

## save RData object
save(eval.o, file=args[7])

