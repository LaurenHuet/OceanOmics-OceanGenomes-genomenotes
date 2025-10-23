module load nextflow/24.04.3

nextflow run main.nf \
   -profile singularity \
   --input assets/samplesheet.csv \
   --outdir results \
   --taxdump /data/databases/taxdump \
   --blastp /data/databases/diamond \
   --blastn /data/databases/blastn