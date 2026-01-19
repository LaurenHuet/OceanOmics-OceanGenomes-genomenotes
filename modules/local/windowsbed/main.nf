process WINDOWSBED {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(fasta)
    val get_sizes

    output:
    tuple val(meta), path("*.{fa,fasta}")     , emit: fa, optional: true
    tuple val(meta), path("*.sizes")          , emit: sizes, optional: true
    tuple val(meta), path("*.fai")            , emit: fai, optional: true
    tuple val(meta), path("*.gzi")            , emit: gzi, optional: true
    tuple val(meta), path("windows.1k.bed")   , emit: windows_bed
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def get_sizes_command = get_sizes ? "cut -f 1,2 ${fasta}.fai > ${fasta}.sizes" : ''

    """
    # Build FASTA index (creates ${fasta}.fai and possibly ${fasta}.gzi)
    samtools \\
        faidx \\
        $fasta \\
        $args

    ${get_sizes_command}

    # Use the newly generated index
    FAI_FILE="${fasta}.fai"

    # Make 1kb windows BED (no overlap)
    cut -f1,2 "\$FAI_FILE" \\
      | awk -v W=1000 'BEGIN{OFS="\\t"} {for(i=0;i<\$2;i+=W){end=i+W; if(end>\$2) end=\$2; print \$1,i,end}}' \\
      > windows.1k.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d;s/.* //')
    END_VERSIONS
    """

    stub:
    def get_sizes_command = get_sizes ? "touch ${fasta}.sizes" : ''
    """
    touch ${fasta}.fai
    if [[ "${fasta.extension}" == "gz" ]]; then
        touch ${fasta}.gzi
    fi

    ${get_sizes_command}
    touch windows.1k.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d;s/.* //')
    END_VERSIONS
    """
}