// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


include { YAML } from '../../../modules/local/yaml/main'
include { CAT_HIC } from '../../../modules/local/cat_hic/main'
include { CAT_HIFI } from '../../../modules/local/cat_hifi/main'
include { FASTP } from '../../../modules/nf-core/fastp/main'
include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HIFI } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HIC } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_HIC } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_HIFI } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { BLOBTK_DEPTH } from '../../../modules/nf-core/blobtk/depth/main'
include { SEQKIT_SLIDING } from '../../../modules/nf-core/seqkit/sliding/main'
include { SEQKIT_SPLIT2 } from '../../../modules/nf-core/seqkit/split2/main'
include { BLAST_BLASTN } from '../../../modules/nf-core/blast/blastn/main'
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
    samplesheet_ch
        .map { meta, hifi_path, hic_path, assembly_path, busco_path ->
            // Extract the required values from the meta map
            def tolid = meta.tolid
            def bioproject_id = meta.bioproject_id
            def taxid = meta.taxid
            def species = meta.species
            
            return [ meta, tolid, bioproject_id, taxid, species ]
        }
        .set { ch_yaml_input }
    
    YAML(ch_yaml_input)
    ch_versions = ch_versions.mix(YAML.out.versions)

    //
    // MODULE: Concatenate HiC reads
    //
    ch_hic = samplesheet_ch
        .map { meta, hifi_path, hic_path, assembly_path, busco_path ->
            def hic_files = file("${hic_path}/*.{fastq,fq,fastq.gz,fq.gz}").findAll { it.exists() } 
            if (hic_files.isEmpty()) {
                error "No HiC files found in ${hic_path}"
            }
            return [ meta, hic_files ]
        }

    CAT_HIC(ch_hic)
    ch_versions = ch_versions.mix(CAT_HIC.out.versions.first())

        //
    // MODULE: Quality filter HiC reads with FASTP
    //
    FASTP(
        CAT_HIC.out.cat_files,
        [], // adapter_fasta - empty for auto-detection
        [], // save_trimmed_fail - empty
        [], // save_merged - empty
        [], // skip_fastp - empty
        "hic" // suffix for output files
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())


    //
    // MODULE: Concatenate HiFi reads
    //
    ch_hifi = samplesheet_ch
        .map { meta, hifi_path, hic_path, assembly_path, busco_path ->
            // Get HiFi files from the directory
            def hifi_files = file("${hifi_path}/*.{fastq,fq,fastq.gz,fq.gz}").findAll { it.exists() }
            if (hifi_files.isEmpty()) {
                error "No HiFi files found in ${hifi_path}"
            }
            return [ meta, hifi_files ]
        }

    CAT_HIFI(ch_hifi)
    ch_versions = ch_versions.mix(CAT_HIFI.out.versions.first())

        //
    // Prepare assembly channel
    //
    ch_assembly = samplesheet_ch
        .map { meta, hifi_path, hic_path, assembly_path, busco_path ->
            // Get the actual assembly file from the directory
            def assembly_files = file("${assembly_path}/*.{fa,fasta,fna}").findAll { it.exists() }
            if (assembly_files.isEmpty()) {
                error "No assembly files found in ${assembly_path}"
            }
            // Take the first assembly file found
            return [ meta, assembly_files[0] ]
        }

    //
    // MODULE: Index the assembly with MINIMAP2
    //
    MINIMAP2_INDEX(ch_assembly)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())

    //
    // MODULE: Align HiFi reads to assembly
    //
    MINIMAP2_ALIGN_HIFI(
        CAT_HIFI.out.hifi_cat,
        ch_assembly,
        MINIMAP2_INDEX.out.index,
        [], // bed - empty
        [], // cigar_paf - empty
        [], // cigar_bam - empty
        "hifi" // suffix
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HIFI.out.versions.first())

    //
    // MODULE: Align HiC reads to assembly
    //
    MINIMAP2_ALIGN_HIC(
        FASTP.out.fastp_hic,
        ch_assembly,
        MINIMAP2_INDEX.out.index,
        [], // bed - empty
        [], // cigar_paf - empty
        [], // cigar_bam - empty
        "hic" // suffix
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HIC.out.versions.first())

    //
    // MODULE: Sort HiFi BAM files
    //
    SAMTOOLS_SORT_HIFI(
        MINIMAP2_ALIGN_HIFI.out.bam,
        ch_assembly.map { meta, fasta -> fasta }.first(),
        "bai"
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_HIFI.out.versions.first())

    //
    // MODULE: Sort HiC BAM files
    //
    SAMTOOLS_SORT_HIC(
        MINIMAP2_ALIGN_HIC.out.bam,
        ch_assembly.map { meta, fasta -> fasta }.first(),
        "bai"
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_HIC.out.versions.first())

//
// MODULE: Merge HiFi and HiC BAM files for combined coverage analysis
//
//
// MODULE: Merge HiFi and HiC BAM files for combined coverage analysis
//
ch_bams_to_merge = SAMTOOLS_SORT_HIFI.out.bam
    .join(SAMTOOLS_SORT_HIC.out.bam)
    .map { meta, hifi_bam, hic_bam ->
        return [ meta, [hifi_bam, hic_bam] ]
    }

SAMTOOLS_MERGE(
    ch_bams_to_merge,
    ch_assembly.map { meta, fasta -> fasta }.first()
)

//
// MODULE: Index the merged BAM file
//
SAMTOOLS_INDEX(SAMTOOLS_MERGE.out.bam)
ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

//
// MODULE: Calculate depth from merged BAM file using BLOBTK_DEPTH
// Uses 1000 bp windows as per Sanger approach - combines HiFi + HiC coverage
//
ch_merged_bam_with_index = SAMTOOLS_MERGE.out.bam
    .join(SAMTOOLS_INDEX.out.bai)

BLOBTK_DEPTH(ch_merged_bam_with_index)
ch_versions = ch_versions.mix(BLOBTK_DEPTH.out.versions.first())

    //
    // MODULE: SEQKIT_SLIDING - Create sliding windows from assembly for contamination screening
    // Using Sanger-ToL BlobToolKit parameters: 100kb chunks, no overlap
    //
    SEQKIT_SLIDING(ch_assembly)
    ch_versions = ch_versions.mix(SEQKIT_SLIDING.out.versions.first())

    //
    // MODULE: SEQKIT_SPLIT2 - Split sliding window output into chunks for parallel BLAST
    // Max 10 chunks per scaffold, min 1000bp length as per Sanger approach
    //
    SEQKIT_SPLIT2(SEQKIT_SLIDING.out.sliding)
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions.first())

    //
    // Prepare channels for BLAST - flatten the split files and maintain meta.id mapping AND taxid
    //
    ch_blast_input = SEQKIT_SPLIT2.out.split_fastx
        .transpose() // This separates each file in the collection into individual channel items
        .map { meta, split_file ->
            // Create unique meta for each split file to ensure separate jobs
            // BUT preserve the original taxid and other important metadata
            def split_meta = [
                id: "${meta.id}_${split_file.baseName}",
                original_id: meta.id,
                taxid: meta.taxid,  // Preserve the taxid from original meta
                species: meta.species ?: "", // Preserve species if available
                tolid: meta.tolid ?: "",
                bioproject_id: meta.bioproject_id ?: ""
            ]
            return [ split_meta, split_file ]
        }

    //
    // MODULE: BLAST_BLASTN - Run BLAST on each chunk with contamination screening
    // Each split file will run as a separate job with taxid exclusion for contamination screening
    //
    BLAST_BLASTN(
        ch_blast_input,
        params.blastn, // Pass database path directly
        [], // taxidlist - empty if not used
        ch_blast_input.map { meta, file -> meta.taxid }, // Use taxid from meta for exclusion
        true  // negative_tax - set to true to exclude the taxid (contamination screening)
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    //
    // // Collect BLAST results per sample
    // // Group by original_id to combine all BLAST results for each sample
    // //
    // ch_blast_results_grouped = BLAST_BLASTN.out.txt
    //     .map { meta, blast_file ->
    //         // Create a new meta with the original_id as the main id for grouping
    //         def group_meta = [
    //             id: meta.original_id,
    //             taxid: meta.taxid,
    //             species: meta.species,
    //             tolid: meta.tolid,
    //             bioproject_id: meta.bioproject_id
    //         ]
    //         return [ group_meta, blast_file ]
    //     }
    //     .groupTuple(by: 0) // Group by meta (which now contains the original_id)
    //     .map { meta, blast_files ->
    //         // blast_files is now a list of all BLAST result files for this sample
    //         return [ meta, blast_files ]
    //     }

    emit:
    yaml_files = YAML.out.yml                           // channel: [ meta, yml_file ]
    hic_reads = CAT_HIC.out.cat_files                   // channel: [ meta, [R1, R2] ]
    hic_reads_filtered = FASTP.out.fastp_hic                // channel: [ meta, [R1_filtered, R2_filtered] ]
    hifi_reads = CAT_HIFI.out.hifi_cat                  // channel: [ meta, hifi_fastq ]
    assembly = ch_assembly                              // channel: [ meta, assembly_fasta ]
    minimap2_index = MINIMAP2_INDEX.out.index          // channel: [ meta, index ]
    hifi_bam_sorted = SAMTOOLS_SORT_HIFI.out.bam       // channel: [ meta, sorted_bam ]
    hic_bam_sorted = SAMTOOLS_SORT_HIC.out.bam         // channel: [ meta, sorted_bam ]
    depth_coverage = BLOBTK_DEPTH.out.coverage_bed         // channel: [ meta, coverage_file ]
    sliding_windows = SEQKIT_SLIDING.out.sliding       // channel: [ meta, sliding_fasta ]
    split_chunks = SEQKIT_SPLIT2.out.split_fastx       // channel: [ meta, [split_files] ]
    //blast_results = ch_blast_results_grouped           // channel: [ meta, [blast_result_files] ]
    samplesheet = samplesheet_ch                       // channel: [ meta, hifi_path, hic_path, assembly_path, busco_path ]
    versions = ch_versions  

}