To streamline the process of evaluating multiple methods/VCFs, we use [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html).
The *sveval* docker image that we provide through [DockerHub](https://hub.docker.com/r/jmonlong/sveval/) contains all dependencies necessary to run the snakemake workflow (R+sveval, snakemake, bcftools, bgzip/tabix).

It's a snakemake workflow so it should be as simple as naming the input correctly and running the `snakemake` command.
The following instructions explain how to name the input file and set up the `config.yaml` file.

The VCFs to analyze are placed in a `vcf` folder. 
They must be named following the template: `{exp}-{method}-{sample}.vcf.gz` where:

- `exp` is the experiment name, usually the name of the catalog or gold-standard dataset. For example `hgsvc`.
- `method` is the label of the method used. E.g. `vg` or `bayestyper`.
- `sample` is the sample name and should also match the sample name in the calls and truth-set.

In addition to the calls, the `vcf` folder must also contains the truth-set, named: `{exp}-truth-baseline.vcf.gz`.
This truth-set must contains genotypes for the samples analyzed.

For example, we could place in the `vcf` folder the following files:

- `vcf/hgsvc-vg-HG00514.vcf.gz`
- `vcf/hgsvc-vg-NA19240.vcf.gz`
- `vcf/hgsvc-bayestyper-NA19240.vcf.gz`
- `vcf/hgsvc-bayestyper-HG00514.vcf.gz`
- `vcf/hgsvc-truth-baseline.vcf.gz`

With these VCF files we can evaluate both methods in two samples using the specified truth-set using:

```
snakemake --configfile config.yaml --cores 8
```

This command reads information from a `config.yaml` file and run the commands in parallel using 8 cores.

Of note, any parameter in the config file can be overwritten in the command line  `--config`.

### Reference genome

The VCF are normalized before evaluation. 
This step requires a FASTA file of the reference genome. 
The path to the reference file must be specified using the config parameter *ref_fa* (see [`config.yaml` file](config.yaml)).

### Non-repeat regions

We usually evaluate the SVs both across the whole-genome and in non-repeat regions only.
The non-repeat regions are defined in a BED file whose path must be specified in the config parameter *nonrep_bed* (see [`config.yaml` file](config.yaml)).

For GRCh38 we use this file [hg38_non_repeats.bed.gz](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/sveval/hg38_non_repeats.bed.gz).

### Output files

The output of this pipeline includes:

- a PDF with bar graph and precision-recall curves in *out_pdf* (defined in  [`config.yaml` file](config.yaml))
- a merged TSV with number of FP/FN/TP and F1, precision, recall scores, in the `tsv` folder.

### Dependencies

In addition to [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), the following must be installed:

- R with :
    - the **sveval** package
    - the **ggplot2** package
    - the **dplyr** package
- [bcftools](https://samtools.github.io/bcftools/bcftools.html)
