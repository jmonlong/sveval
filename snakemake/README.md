To streamline the process of evaluating multiple methods/VCFs, we use [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html).
The *sveval* docker image that we provide in [quay.io/jmonlong/sveval](https://quay.io/repository/jmonlong/sveval) contains all dependencies necessary to run the snakemake workflow (R+sveval, snakemake, bcftools, bgzip/tabix).
For example, use `quay.io/jmonlong/sveval:v2.2.0` (see [example below](#start-the-docker-container)).

It's a snakemake workflow so it should be as simple as naming the input correctly and running the `snakemake` command.

The following instructions explain how to name the input files and set up the `config.yaml` file.

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

For GRCh38 we use this file: [hg38_non_repeats.bed.gz](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/sveval/hg38_non_repeats.bed.gz).

### Simple repeat track

To help match SVs within simple repeats, we can provide a BED file of the simple repeat annotation.
sveval will allow more wiggle room, i.e. for a SV to "move" along a simple repeat segment.
It helps matching similar simple repeat variants that are just positioned differently.
For example, a similar deletion might be called at the beginning of the annotated repeat but positioned at the end in the truthset.
And because simple repeats are not perfect and the SVs exactly the same, left-aligning might not handle these cases.

In the config file, the BED file with the simple repeat annotation is specified with the *simprep_bed* parameter (see [`config.yaml` file](config.yaml)).
Comment that line out to disable this feature (although it is highly recommended to use it).

Any BED file would work, but we've also prepared some for GRCh38 and GRCh37: [../docs/simpleRepeat_GRCh38.bed.gz](../docs/simpleRepeat_GRCh38.bed.gz) and [../docs/simpleRepeat_GRCh37.bed.gz](../docs/simpleRepeat_GRCh37.bed.gz).
To match the config in [`config.yaml` file](config.yaml)/[`config.example.yaml`](config.example.yaml), those files would be placed in a *bed* folder.
Of note, they've been prepared with by downloading the TRF annotation from the UCSC Genome Browser and trimming/merging information (see [this notebook](../docs/prepare-simple-repeat-track-from-ucsc.md)).


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

As mentioned above, a docker container with all the dependencies is available at [quay.io/jmonlong/sveval](https://quay.io/repository/jmonlong/sveval) (e.g. `quay.io/jmonlong/sveval:v2.2.0`). 
See [example below](#start-the-docker-container).


### Example: evaluating SVs on GRCh37 agains the GIAB SV truthset

Here are the commands one could use to evaluate a VCF file against the GIAB truthset from scratch.

#### Clone the repo and move to the snakemake pipeline directory

```sh
git clone https://github.com/jmonlong/sveval.git
cd sveval/snakemake
```

#### Download (or copy) the reference FASTA

```sh
wget -O hs37d5.fa.gz ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
```

#### Place VCFs in a 'vcf' directory

```sh
mkdir vcf
```

Download VCF with the truthset from GIAB:

```sh
wget -O vcf/giab-truth-baseline.vcf.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
```

Place the VCFs with the calls in the vcf folder.
In this example, SV calls using vg on HG002 and GRCh37:

```sh
wget -O vcf/giab-vg-HG002.vcf.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/vgsv2019/vcfs/giab5-vg-HG002.vcf.gz
```

#### Download the BED file defining confident regions for this truth set

```sh
mkdir bed
wget -O bed/conf.bed ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
```

#### Download the BED file defining simple repeats

```sh
wget -O bed/simpleRepeat_GRCh37.bed.gz https://raw.githubusercontent.com/jmonlong/sveval/master/docs/simpleRepeat_GRCh37.bed.gz
```

#### Edit config.yaml to match the files/experiment

In this example, we evaluate only the 'vg' method, on the 'HG002' sample, and for the 'giab' experiment/truthset.
We also have a different reference fasta, regions BED, and simple repeat BED.
See [`config.example.yaml`](config.example.yaml).

#### Start the docker container

- `-v`/`-w` to bind the current directory to a /app working directory in the container
- `-u` to make sure files are created with permissions matching the user (and not created by 'root') 

```sh
docker run -it --rm -v `pwd`:/app -w /app -u `id -u $USER` quay.io/jmonlong/sveval:v2.2.0
```

#### Run the snakemake pipeline within the container

```sh
snakemake --configfile config.example.yaml --cores 4
```
