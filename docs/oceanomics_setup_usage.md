## Overview of workflow

Before running the Nextflow pipeline itself, three preparatory steps are required:

1. Configure the pipeline using a single configuration file
2. Generate a samplesheet from the OceanOmics PostgreSQL database
3. Stage required input data (HiFi, Hi-C, assemblies, BUSCO results)

The generated samplesheet is the **single source of truth** for all downstream steps.

---

## Pre-run setup

### 1. Configure the pipeline

All preparatory steps are controlled by a bash-compatible configuration file:

```
genomenotes_pipeline.conf
```

This file defines which samples to include, where data should be staged, and where database credentials are stored.

Minimum required fields:

```
# Path to Postgres credentials file
POSTGRES_CFG=~/postgresql_details/oceanomics.cfg

# OceanOmics sample identifiers to include
OG_IDS="OG38"

# Base directory where all data will be staged
# {user} is replaced automatically at runtime
STAGING_BASE_DIR=/scratch/pawsey0964/{user}/

# Where the generated samplesheet will be written
SAMPLESHEET_OUTPUT_DIR=/scratch/pawsey0964/{user}/
SAMPLESHEET_FILENAME_PREFIX=samplesheet

# Path to Samplesheet to be used by staging scripts
SAMPLESHEET=/scratch/pawsey0964/{user}/
```

PostgreSQL credentials are **not** stored in this file. Instead, `POSTGRES_CFG` must point to an INI-style file with the following format:

```
[postgres]
dbname = oceanomics
user = <username>
password = <password>
host = <hostname>
port = 5432
```

---

### 2. Generate the samplesheet

The samplesheet is generated directly from the OceanOmics PostgreSQL database and contains:

- Sample identifiers
- Paths to staged data directories
- Assembly versioning information
- Taxonomic metadata

To generate the samplesheet:

```
python scripts/create_samplesheet_from_config.py ../genomenotes_pipeline.conf
```

This will create a dated CSV file (e.g. `samplesheet_20260107.csv`) in the directory specified by `SAMPLESHEET_OUTPUT_DIR`.

If paths need to be changed, update the configuration file and **regenerate the samplesheet**.

---

### 3. Stage input data

Input data are staged using helper scripts, all driven by the same configuration file and samplesheet.

Each script reads the samplesheet and copies data from the appropriate rclone remotes into the per-sample directories defined in the samplesheet.

Stage HiFi reads:

```
bash scripts/02_get_hifi_from_config.sh ../genomenotes_pipeline.conf
```

Stage Hi-C reads:

```
bash scripts/01_get_hic_from_config.sh ../genomenotes_pipeline.conf
```

Stage assemblies:

```
bash scripts/03_get_assembly_from_config.sh ../genomenotes_pipeline.conf
```

Stage BUSCO results:

```
bash scripts/04_get_busco_from_config.sh ../genomenotes_pipeline.conf
```

These scripts are safe to re-run. Existing files will be skipped by `rclone`.

After staging, each sample directory should contain:

```
OG38/
├── hifi/
├── hic/
├── assembly/
└── busco/
```