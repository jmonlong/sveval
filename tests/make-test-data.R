library(VariantAnnotation)

##
## add some genotype quality to the multi-samp VCF
##
vcf = readVcf('truth.refalt.vcf')

## GQ matrix
gt = geno(vcf)$GT
gq = matrix(round(runif(length(gt), 0, 30)), nrow(gt), ncol(gt))
rownames(gq) = rownames(gt)
colnames(gq) = colnames(gt)
geno(vcf)$GQ = gq

## header
hh = geno(header(vcf))
hh = rbind(hh, S4Vectors::DataFrame(Number='1', Type='Integer', Description='Genotype Quality'))
rownames(hh) = c('GT', 'GQ')
geno(header(vcf)) = hh

writeVcf(vcf, 'ex.gq.vcf')
