process LIFT_COORDS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(blast_files)

    output:
    tuple val(meta), path("*.lifted.tsv"), emit: tsv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create directory and copy files
    mkdir -p cat_blast
    cp ${blast_files.join(' ')} cat_blast/
    
    # Concatenate all BLAST result files
    cat cat_blast/*.txt > combined_blastn_results.tsv

    # Run the lift_coords.py script
    lift_coords.py \\
        combined_blastn_results.tsv \\
        ${prefix}.lifted.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lifted.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}