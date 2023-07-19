## Documentation

## Manual

- [Manual (PDF version)](sveval-manual.pdf)

## Preparing the simple repeat track

To help match SVs within simple repeats, we can provide a BED file of the simple repeat annotation.
sveval will allow more wiggle room, i.e. for a SV to "move" along a simple repeat segment.
It helps matching similar simple repeat variants that are just positioned differently.
For example, a similar deletion might be called at the beginning of the annotated repeat but positioned at the end in the truthset.
And because simple repeats are not perfect and the SVs exactly the same, left-aligning might not handle these cases.

We've prepared simple repeats regions for GRCh38 and GRCh37: [simpleRepeat_GRCh38.bed.gz](simpleRepeat_GRCh38.bed.gz) and [simpleRepeat_GRCh37.bed.gz](simpleRepeat_GRCh37.bed.gz).
More details on this in [*prepare-simple-repeat-track-from-ucsc.md* markdown report](prepare-simple-repeat-track-from-ucsc.md).
