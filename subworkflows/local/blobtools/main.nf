// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


include {YMAL} from '../../../modules/local/yaml/main'
include { FASTP } from '../../../modules/nf-core/fastp/main'
include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { BLOBTK_DEPTH } from '../../../modules/nf-core/blobtk/depth/main'
include { SEQKIT_SLIDING } from '../../../modules/nf-core/seqkit/sliding/main'
include { SEQKIT_SPLIT2 } from '../../../modules/nf-core/seqkit/split2/main'
include { BLAST_BLASTN } from '../../../modules/nf-core/blast/blastn/main'
include { BLOBTOOLS_ADD } from '../../../modules/local/blobtools/add/main'
include { BLOBTOOLS_CREATE } from '../../../modules/local/blobtools/create/main'

workflow BLOBTOOLS {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
