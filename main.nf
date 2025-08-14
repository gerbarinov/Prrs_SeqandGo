#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMMOMATIC } from './modules.nf'
include { BWA_INDEX } from './modules.nf'
include { BWA_MEM } from './modules.nf'
include { INDEX_BAM } from './modules.nf'
include { COLLECT_MAPPING_STATS } from './modules.nf'
include { CREATE_CONSENSUS } from './modules.nf'
include { NEXTSTRAIN_PHYLOGENY } from './modules.nf'

// Define input channels
read_pairs_ch = Channel.fromFilePairs("/home/ubuntu/data/PRRS/*_L00_{R1,R2}.fq.gz", checkIfExists: true)
references_ch = Channel.fromList(params.references.collect { ref ->
    def refFile = file(ref)
    [refFile.simpleName, refFile]
})

workflow {
    // Trimming
    TRIMMOMATIC(read_pairs_ch, params.adapters, params.primers)

    // Index references
    BWA_INDEX(references_ch)

    // Combine trimmed reads with indexed references
    trimmed_with_refs = TRIMMOMATIC.out.trimmed_reads
        .combine(BWA_INDEX.out.index)
        .map { sample_id, reads, ref_name, ref_file, index_files ->
            [sample_id, reads, ref_name, ref_file, index_files.flatten()]
        }

    // Mapping to both references
    BWA_MEM(trimmed_with_refs)

    // Index BAM files
    INDEX_BAM(BWA_MEM.out.bam)

    // Collect mapping statistics
    COLLECT_MAPPING_STATS(INDEX_BAM.out.indexed_bam)

    // Combine indexed BAM files with references for consensus creation
    bam_ref_combined = INDEX_BAM.out.indexed_bam
        .combine(references_ch)
        .map { sample_id, ref_name_bam, bam, bai, ref_name_fasta, fasta ->
            [sample_id, ref_name_fasta, bam, bai, fasta]
        }

    // Create consensus sequence
    CREATE_CONSENSUS(bam_ref_combined)

    // Collect and group consensus sequences
    CREATE_CONSENSUS.out.consensus
        .filter { it[2].size() > 0 }
        .set { non_empty_consensus }

    // Collect all non-empty consensus files
    non_empty_consensus
        .map { it[2] }
        .collect()
        .set { collected_consensus }

    // Run Nextstrain phylogenetic analysis
    NEXTSTRAIN_PHYLOGENY(
        collected_consensus,
        file(params.nextstrain_config)
    )

    // View results
    INDEX_BAM.out.indexed_bam.view { "Indexed BAM: $it" }
    COLLECT_MAPPING_STATS.out.mapping_stats.view { "Mapping stats: $it" }
    non_empty_consensus.view { "Non-empty consensus file: ${it[2]}" }
    NEXTSTRAIN_PHYLOGENY.out.nextstrain_results.view { "Nextstrain results: $it" }
}