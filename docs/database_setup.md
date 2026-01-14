## Overview of database setups

You will need access to the following databases 

### 1 NCBI taxdump database 

Create the database directory, retreive and decompress the NCBI taxonomy:

```
DATE=2026_01_07
TAXDUMP=/path/to/databases/taxdump_${DATE}
TAXDUMP_TAR=/path/to/databases/taxdump_${DATE}.tar.gz
mkdir -p "$TAXDUMP"
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz -o $TAXDUMP_TAR
tar -xzf $TAXDUMP_TAR -C "$TAXDUMP"
```

### 2 NCBI core nucleotide BLAST database

This can be retrevied from /scratch/references/ on pawsey, check for the latest verion of the database and copy into the NT directory. 

Create the database directory and move into the directory:

```
DATE=2026_01
NT=/scratch/pawsey0964/lhuet/PIPELINE_DEV/GEOMENOTES/core_nt_${DATE}
NT_TAR=/scratch/pawsey0964/lhuet/PIPELINE_DEV/GEOMENOTES/core_nt_${DATE}.tar.gz
mkdir -p $NT
cd /scratch/references/blastdb_update/blast-2026-01-05 
## make sure you get the latest database
rclone copy core_nt.* $NT -v
```