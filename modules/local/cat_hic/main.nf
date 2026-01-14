process CAT_HIC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(hic_files)

    output:
    tuple val(meta), path("cat_files/*_R*.fastq.gz"), emit: cat_files
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    # If multiple lane/files exist, concatenate them; otherwise just rename consistently
    if [ "\$(ls *R1*fastq.gz 2>/dev/null | wc -l)" -gt 1 ]; then
        cat \\
            $args \\
            *R1*fastq.gz \\
            > ${prefix}_R1.fastq.gz

        cat \\
            $args \\
            *R2*fastq.gz \\
            > ${prefix}_R2.fastq.gz
    else
        mv *R1*fastq.gz ${prefix}_R1.fastq.gz
        mv *R2*fastq.gz ${prefix}_R2.fastq.gz
    fi

    mkdir cat_files
    mv ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz cat_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
