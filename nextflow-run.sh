#!/bin/bash

module load nextflow/24.10.0

export NXF_HOME=/scratch/pawsey0964/lhuet/.nextflow
export NXF_SINGULARITY_CACHEDIR="/scratch/pawsey0964/lhuet/singularity/.nextflow_singularity"
mkdir -p $NXF_HOME/plugins
mkdir -p $NXF_SINGULARITY_CACHEDIR

echo "Set NXF_HOME to: $NXF_HOME"
echo "Set NXF_SINGULARITY_CACHEDIR to: $NXF_SINGULARITY_CACHEDIR"
echo "Created plugins directory: $NXF_HOME/plugins"

nextflow run main.nf \
   -profile singularity \
   --input assets/samplesheet.csv \
   --outdir /scratch/pawsey0964/lhuet/PIPELINE_DEV/GEOMENOTES \
   --taxdump /scratch/pawsey0964/lhuet/PIPELINE_DEV/GEOMENOTES/taxdump_2026_01_07 \
   --blastp /scratch/pawsey0964/lhuet/databases/RefSeq_protein_25_10_01 \
   --blastn /scratch/pawsey0964/lhuet/PIPELINE_DEV/GEOMENOTES/blastdb_14_01 \
   -c pawsey_profile.config -resume \
   -c local_blobtoolkit.config \
   -with-report usage_report_$(date +%Y%m%d_%H%M%S).html \
   --binddir /scratch 