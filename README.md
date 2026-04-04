# LaurenHuet/OceanOmics-OceanGenomes-genomenotes

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/LaurenHuet/OceanOmics-OceanGenomes-genomenotes)
[![GitHub Actions CI Status](https://github.com/LaurenHuet/OceanOmics-OceanGenomes-genomenotes/actions/workflows/nf-test.yml/badge.svg)](https://github.com/LaurenHuet/OceanOmics-OceanGenomes-genomenotes/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/LaurenHuet/OceanOmics-OceanGenomes-genomenotes/actions/workflows/linting.yml/badge.svg)](https://github.com/LaurenHuet/OceanOmics-OceanGenomes-genomenotes/actions/workflows/linting.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.4.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064)](https://github.com/nf-core/tools/releases/tag/3.4.1)

---

## Introduction

**OceanOmics-OceanGenomes-genomenotes** is an nf-core–style pipeline for generating genome notes for the Ocean Genomes project, following the principles and structure of the Sanger Tree of Life genome notes. The pipeline integrates sample metadata from the OceanOmics database with staged sequencing data and assemblies to produce standardised, reproducible genome notes.

---

## OceanOmics internal useage

To set up the pipeline follow the steps in

- [Oceanomics usage](docs/oceanomics_setup_usage.md)
  - Contains information on how to use scripts to set up samplesheet from database and how to stage all data using a central configuration file.
---

## Running the pipeline

The pipeline requires a samplesheet with the following format

`samplesheet.csv`:

```csv
sample,hifi_reads,hic_reads,assembly,busco_genes,bioproject_id,version,date,tolid,taxid,species
sample1,/path/to/hifi_reads_dir,/path/to/hic_reads_dir,/path/to/assembly_dir,/path/to/busco_genes_dir/bioproject_id,version,date,tolid,taxid,species
sample2,/path/to/hifi_reads_dir2,/path/to/hic_reads_dir2,/path/to/assembly_dir2,/path/to/busco_genes_dir2/bioproject_id,version,date,tolid,taxid,species
```
Each row represents a sample with the following columns:
- `sample`: Sample identifier (must be unique)
- `hifi_reads`: Path to directory containing HiFi reads
- `hic_reads`: Path to directory containing HiC reads  
- `assembly`: Path to directory containing assembly files
- `busco_genes`: Path to directory containing BUSCO genes
- `bioproject_id`: From the uploaded genome
- `version`: Oceanomics assembly version
- `date`: Date of assembly
- `tolid`: Tree of Life ID for sample
- `taxid`: NCBI taxon ID for sample
- `species` : Validated Species ID

You will also need to prepare the following databases:
- **Taxdump database**: For taxonomic classification
- **Diamond BLASTP database**: For protein sequence alignment
- **BLASTN database**: For nucleotide sequence alignment


Once the samplesheet has been generated and all data are staged, the pipeline can be run using Nextflow:



```bash
nextflow run LaurenHuet/OceanOmics-OceanGenomes-genomenotes \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --taxdump /path/to/taxdump/database \
   --blastp /path/to/diamond/database \
   --blastn /path/to/blastn/database
```

### Required Parameters

- `--input`: Path to comma-separated file containing information about the samples
- `--outdir`: The output directory where the results will be saved
- `--taxdump`: Path to the taxdump database directory for taxonomic classification
- `--blastp`: Path to the Diamond BLASTP database directory for protein sequence alignment
- `--blastn`: Path to the BLASTN database directory for nucleotide sequence alignment

### Example Usage

```bash
# Run with Singularity
nextflow run LaurenHuet/OceanOmics-OceanGenomes-genomenotes \
   -profile singularity \
   --input assets/samplesheet.csv \
   --outdir results \
   --taxdump /data/databases/taxdump \
   --blastp /data/databases/diamond \
   --blastn /data/databases/blastn
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

LaurenHuet/OceanOmics-OceanGenomes-genomenotes was originally written by LaurenHuet.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use LaurenHuet/OceanOmics-OceanGenomes-genomenotes for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
a