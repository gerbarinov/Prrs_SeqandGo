#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMMOMATIC } from './modules.nf'
include { SPADES_ASSEMBLY_RESCUE } from './modules.nf'
include { EXTRACT_ORF_FROM_ASSEMBLY } from './modules.nf'
include { BWA_INDEX } from './modules.nf'
include { BWA_MEM } from './modules.nf'
include { INDEX_BAM } from './modules.nf'
include { DOWNSAMPLE_BAM } from './modules.nf'
include { COLLECT_MAPPING_STATS } from './modules.nf'
include { CREATE_CONSENSUS } from './modules.nf'
include { NEXTSTRAIN_PHYLOGENY } from './modules.nf'

// Define input channels - FIXED pattern to match config
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

// Create references channel for prrsV1 and prrsV2 only
references_ch = Channel.fromList(params.references.collect { ref ->
    def refFile = file(ref)
    [refFile.simpleName, refFile]
})

workflow {
    // Trimming
    TRIMMOMATIC(read_pairs_ch, params.adapters, params.primers)

    // Index the two main references (prrsV1 and prrsV2)
    BWA_INDEX(references_ch)

    // Combine trimmed reads with indexed references
    trimmed_with_refs = TRIMMOMATIC.out.trimmed_reads
        .combine(BWA_INDEX.out.index)
        .map { sample_id, reads, ref_name, ref_file, index_files ->
            [sample_id, reads, ref_name, ref_file, index_files.flatten()]
        }

    // Mapping to both references (prrsV1 and prrsV2)
    BWA_MEM(trimmed_with_refs)

    // Index BAM files
    INDEX_BAM(BWA_MEM.out.bam)

    // Downsample BAM files to 100,000 reads
    DOWNSAMPLE_BAM(INDEX_BAM.out.indexed_bam)

    // Collect mapping statistics on downsampled BAM files
    COLLECT_MAPPING_STATS(DOWNSAMPLE_BAM.out.downsampled_bam)

    // Combine downsampled BAM files with references for consensus creation
    bam_ref_combined = DOWNSAMPLE_BAM.out.downsampled_bam
        .combine(references_ch)
        .map { sample_id, ref_name_bam, bam, bai, ref_name_fasta, fasta ->
            [sample_id, ref_name_fasta, bam, bai, fasta]
        }

    // Create consensus sequence for ORF5 and ORF7 regions
    CREATE_CONSENSUS(bam_ref_combined)

    // Check which samples need rescue assembly
    // Group consensus results by sample to check if any reference has good consensus
    consensus_by_sample = CREATE_CONSENSUS.out.consensus
        .map { sample_id, ref_name, orf5_consensus, orf7_consensus ->
            [sample_id, ref_name, orf5_consensus.size(), orf7_consensus.size()]
        }
        .groupTuple(by: 0)
        .map { sample_id, ref_names, orf5_sizes, orf7_sizes ->
            // Check if any reference has at least one ORF > 2 bytes
            def has_good_consensus = false
            for (int i = 0; i < ref_names.size(); i++) {
                if (orf5_sizes[i] > 2 || orf7_sizes[i] > 2) {
                    has_good_consensus = true
                    break
                }
            }
            [sample_id, has_good_consensus]
        }

    // Get samples that need rescue assembly (no good consensus for any reference)
    samples_needing_rescue = consensus_by_sample
        .filter { sample_id, has_good_consensus -> !has_good_consensus }
        .map { sample_id, has_good_consensus -> sample_id }

    // Get trimmed reads for samples that need rescue assembly
    rescue_reads = samples_needing_rescue
        .combine(TRIMMOMATIC.out.trimmed_reads, by: 0)

    // Run rescue SPAdes assembly for samples with no good consensus (once per sample, no reference template)
    SPADES_ASSEMBLY_RESCUE(rescue_reads)

    // Extract ORF sequences from rescue assemblies for both references
    rescue_assembly_with_refs = SPADES_ASSEMBLY_RESCUE.out.assembly
        .combine(references_ch)
        .map { sample_id, assembly, ref_name, ref_file ->
            [sample_id, assembly, ref_name, ref_file]
        }

    EXTRACT_ORF_FROM_ASSEMBLY(rescue_assembly_with_refs)

    // Combine successful consensus sequences with rescue ORF sequences
    successful_consensus = CREATE_CONSENSUS.out.consensus
        .filter { sample_id, ref_name, orf5_consensus, orf7_consensus ->
            // Keep consensus files that are > 2 bytes
            def orf5_good = orf5_consensus.size() > 2
            def orf7_good = orf7_consensus.size() > 2
            return orf5_good || orf7_good
        }

    // Prepare all sequences for phylogenetic analysis
    all_orf5_sequences = successful_consensus
        .map { sample_id, ref_name, orf5_consensus, orf7_consensus ->
            orf5_consensus
        }
        .filter { it.size() > 2 }
        .mix(
            EXTRACT_ORF_FROM_ASSEMBLY.out.orf5_sequence
                .filter { it.size() > 2 }
        )

    all_orf7_sequences = successful_consensus
        .map { sample_id, ref_name, orf5_consensus, orf7_consensus ->
            orf7_consensus
        }
        .filter { it.size() > 2 }
        .mix(
            EXTRACT_ORF_FROM_ASSEMBLY.out.orf7_sequence
                .filter { it.size() > 2 }
        )

    // Collect all sequences for phylogenetic analysis
    collected_consensus = all_orf5_sequences
        .mix(all_orf7_sequences)
        .collect()

    // Run Nextstrain phylogenetic analysis
    NEXTSTRAIN_PHYLOGENY(
        collected_consensus,
        file(params.nextstrain_config)
    )

    // View results
    DOWNSAMPLE_BAM.out.downsampled_bam.view { "Downsampled BAM: $it" }
    COLLECT_MAPPING_STATS.out.mapping_stats.view { "Mapping stats: $it" }
    successful_consensus.view { "Successful consensus files: ORF5=${it[2]}, ORF7=${it[3]}" }
    samples_needing_rescue.view { "Samples needing rescue assembly: $it" }
    SPADES_ASSEMBLY_RESCUE.out.assembly.view { "Rescue assembly created: $it" }
    EXTRACT_ORF_FROM_ASSEMBLY.out.orf5_sequence.view { "Extracted ORF5 from assembly: $it" }
    EXTRACT_ORF_FROM_ASSEMBLY.out.orf7_sequence.view { "Extracted ORF7 from assembly: $it" }
    NEXTSTRAIN_PHYLOGENY.out.nextstrain_results.view { "Nextstrain results: $it" }
}
