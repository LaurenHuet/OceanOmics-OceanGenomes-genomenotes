process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ 'docker://quay.io/biocontainers/blast:2.17.0--h66d330f_0' }"

    input:
    tuple val(meta) , path(split_fastx)
    val db_prefix_path
    path taxidlist
    val taxids
    val negative_tax

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

script:
def args   = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
def fasta  = split_fastx

def is_compressed = fasta.getExtension() == "gz"
def fasta_name    = is_compressed ? fasta.getBaseName() : fasta.getName()

def negative_flag = negative_tax ? "negative_" : ""

def taxidlist_cmd = (taxidlist) ? "-${negative_flag}taxidlist ${taxidlist}" : ""
def taxids_cmd    = (taxids && taxids != [] && taxids.toString().trim())
                    ? "-${negative_flag}taxids ${taxids}"
                    : ""

if (taxidlist_cmd && taxids_cmd) {
    error "taxidlist and taxids cannot be used at the same time. Choose one."
}

"""


blastn \
    -num_threads ${task.cpus} \
    -db ${db_prefix_path}/core_nt \
    -query ${fasta_name} \
    ${taxidlist_cmd} \
    ${taxids_cmd} \
    ${args} \
    -out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
"""
}
