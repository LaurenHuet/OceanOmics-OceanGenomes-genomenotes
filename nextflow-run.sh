#!/bin/bash

module load nextflow/24.10.0

export NXF_HOME=/scratch/pawsey0964/lhuet/.nextflow
mkdir -p $NXF_HOME/plugins

echo "Set NXF_HOME to: $NXF_HOME"
echo "Created plugins directory: $NXF_HOME/plugins"

nextflow run main.nf \
   -profile singularity \
   --input assets/samplesheet.csv \
   --outdir /scratch/pawsey0964/lhuet/genomenotes \
   --taxdump /scratch/pawsey0964/lhuet/databases/taxdump_25_10 \
   --blastp /scratch/pawsey0964/lhuet/databases/RefSeq_protein_25_10_01 \
   --blastn /scratch/pawsey0964/lhuet/databases/blast_core_nt_2025_10_01 \
   -c pawsey_profile.config -resume \
   --binddir /scratch 