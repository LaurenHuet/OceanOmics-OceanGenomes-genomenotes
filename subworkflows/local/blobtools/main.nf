// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


include { YAML } from '../../../modules/local/yaml/main'
//include { FASTP } from '../../../modules/nf-core/fastp/main'
//include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'
//include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
//include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
//include { BLOBTK_DEPTH } from '../../../modules/nf-core/blobtk/depth/main'
//include { SEQKIT_SLIDING } from '../../../modules/nf-core/seqkit/sliding/main'
//include { SEQKIT_SPLIT2 } from '../../../modules/nf-core/seqkit/split2/main'
//include { BLAST_BLASTN } from '../../../modules/nf-core/blast/blastn/main'
//include { BLOBTOOLS_ADD } from '../../../modules/local/blobtools/add/main'
//include { BLOBTOOLS_CREATE } from '../../../modules/local/blobtools/create/main'

workflow BLOBTOOLS {


    take:
    samplesheet_ch // channel: [ enriched_meta, hifi_path, hic_path, assembly_path, busco_path ]
                   // where enriched_meta contains: id, bioproject_id, version, date, tolid, taxid, species
    
    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Generate YAML configuration files for blobtools
    // Uses metadata from samplesheet: tolid, bioproject_id, taxid, species
    //
    YAML(samplesheet_ch)
    ch_versions = ch_versions.mix(YAML.out.versions)

    // TODO: Add other blobtools processes here
    // The YAML files will be used by subsequent blobtools processes
    // Available channels:
    // - YAML.out.yml: [ meta, yml_file ] - YAML config files
    // - samplesheet_ch: [ meta, hifi_path, hic_path, assembly_path, busco_path ] - original data

    emit:
    yaml_files = YAML.out.yml      // channel: [ meta, yml_file ]
    samplesheet = samplesheet_ch   // channel: [ meta, hifi_path, hic_path, assembly_path, busco_path ]
    versions   = ch_versions       // channel: [ versions.yml ]
}
