process TRIMMOMATIC {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)
    path adapters
    path primers

    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads

    script:
    """
    trimmomatic PE -threads $task.cpus \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_R1_trimmed.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \
        ${sample_id}_R2_trimmed.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${adapters}:2:30:10 \
        ILLUMINACLIP:${primers}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process BWA_INDEX {
    input:
    tuple val(reference_name), path(reference)

    output:
    tuple val(reference_name), path(reference), path("${reference}.*"), emit: index

    script:
    """
    bwa index ${reference}
    """
}

process BWA_MEM {
    input:
    tuple val(sample_id), path(reads), val(reference_name), path(reference), path(index_files)

    output:
    tuple val(sample_id), val(reference_name), path("${sample_id}_${reference_name}.bam"), emit: bam

    script:
    """
    bwa mem -t $task.cpus ${reference} ${reads[0]} ${reads[1]} | \
    samtools view -bS - | \
    samtools sort -o ${sample_id}_${reference_name}.bam -
    """
}

process INDEX_BAM {
    input:
    tuple val(sample_id), val(reference_name), path(bam)

    output:
    tuple val(sample_id), val(reference_name), path("${bam}"), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index ${bam}
    """
}

process COLLECT_MAPPING_STATS {
    input:
    tuple val(sample_id), val(reference_name), path(bam), path(bai)

    output:
    path "${sample_id}_${reference_name}_mapping_stats.txt", emit: mapping_stats

    script:
    """
    samtools flagstat ${bam} > ${sample_id}_${reference_name}_mapping_stats.txt
    samtools idxstats ${bam} >> ${sample_id}_${reference_name}_mapping_stats.txt
    """
}

process CREATE_CONSENSUS {
    conda "bioconda::samtools=1.15 bioconda::bcftools=1.15 bioconda::htslib=1.15"
    
    input:
    tuple val(sample_id), val(ref_name), path(bam), path(bai), path(reference)

    output:
    tuple val(sample_id), val(ref_name), path("${sample_id}_${ref_name}_consensus.fasta"), emit: consensus
    path "${sample_id}_${ref_name}_variants.vcf.gz", emit: variants

    script:
    """
    set -e
    set -x

    REF_NAME=\$(grep '^>' ${reference} | head -n 1 | cut -d ' ' -f 1 | sed 's/>//')
    REGION="\${REF_NAME}:14318-14987"

    # Check if there are any reads in the specified region
    read_count=\$(samtools view ${bam} \$REGION | wc -l)
    echo "Read count in region: \$read_count"

    if [ \$read_count -eq 0 ]; then
        echo "No reads found in the specified region for sample ${sample_id} with reference ${ref_name}. Skipping consensus creation."
        touch ${sample_id}_${ref_name}_consensus.fasta
        touch ${sample_id}_${ref_name}_variants.vcf.gz
    else
        # Generate VCF file for the specific region
        bcftools mpileup -f ${reference} -r \$REGION ${bam} | \
        bcftools call -mv -Oz -o ${sample_id}_${ref_name}_variants.vcf.gz

        # Index the VCF file
        bcftools index ${sample_id}_${ref_name}_variants.vcf.gz

        # Create consensus sequence for the whole genome
        bcftools consensus -f ${reference} -o ${sample_id}_${ref_name}_full_consensus.fasta ${sample_id}_${ref_name}_variants.vcf.gz

        # Extract the region of interest
        samtools faidx ${sample_id}_${ref_name}_full_consensus.fasta \$REGION > ${sample_id}_${ref_name}_consensus.fasta

        # Rename the sequence in the FASTA file
        sed -i "1s/.*/>$sample_id/" ${sample_id}_${ref_name}_consensus.fasta

        # Clean up intermediate files
        rm ${sample_id}_${ref_name}_full_consensus.fasta*
    fi

    echo "Consensus creation completed for sample ${sample_id} with reference ${ref_name}."
    """
}

process NEXTSTRAIN_PHYLOGENY {
    conda "bioconda::nextstrain-cli=3.2.1 bioconda::nextclade=2.14.0"
    
    input:
    path consensus_files
    path nextstrain_config

    output:
    path "results/*", emit: nextstrain_results

    script:
    """
    # Download the PRRS ORF5 database
    nextclade dataset get --name 'community/isuvdl/mazeller/prrsv1/orf5/yimim2025' --output-dir 'datasets/prrsv1'

    # Create a directory for input sequences
    mkdir -p data/sequences

    # Copy consensus files to the input directory
    cp ${consensus_files} data/sequences/

    # Run Nextstrain workflow
    nextstrain build . \
        --configfile ${nextstrain_config} \
        --cores ${task.cpus} \
        --native \
        --workdir ./results

    # Optionally, you can add commands here to process or format the results
    """
}