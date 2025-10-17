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

process SPADES_ASSEMBLY_RESCUE {
    conda "bioconda::spades=3.15.5 conda-forge::python=3.9 conda-forge::setuptools conda-forge::distutils-extra"
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_rescue_assembly/contigs.fasta"), emit: assembly
    path "${sample_id}_rescue_assembly/assembly_graph.fastg", optional: true, emit: graph
    path "${sample_id}_rescue_assembly/spades.log", emit: log

    script:
    """
    set -e
    set -x

    echo "Starting rescue SPAdes assembly for sample ${sample_id} (no reference template)..."
    
    # Check input files
    echo "Checking input files:"
    ls -la ${reads[0]} ${reads[1]}
    
    # Get basic stats about input files
    echo "Input file statistics:"
    r1_reads=\$(zcat ${reads[0]} | wc -l | awk '{print \$1/4}')
    r2_reads=\$(zcat ${reads[1]} | wc -l | awk '{print \$1/4}')
    echo "R1 reads: \$r1_reads"
    echo "R2 reads: \$r2_reads"

    # Run SPAdes with RNA viral mode for better viral genome assembly (no reference template)
    spades.py --rnaviral  \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o ${sample_id}_rescue_assembly \
        --threads ${task.cpus} \
        --memory ${task.memory.toGiga()}

    # Check if assembly was successful
    if [ ! -f "${sample_id}_rescue_assembly/contigs.fasta" ]; then
        echo "ERROR: Rescue SPAdes assembly failed - no contigs.fasta produced"
        exit 1
    fi

    # Report assembly statistics
    echo "Rescue assembly completed for sample ${sample_id}"
    contig_count=\$(grep -c ">" ${sample_id}_rescue_assembly/contigs.fasta || echo "0")
    echo "Number of contigs: \$contig_count"
    """
}

process EXTRACT_ORF_FROM_ASSEMBLY {
    conda "bioconda::blast=2.12.0 bioconda::seqtk=1.3"
    tag "$sample_id-$ref_name"
    
    input:
    tuple val(sample_id), path(assembly), val(ref_name), path(reference)

    output:
    path "${sample_id}_${ref_name}_ORF5_from_assembly.fasta", emit: orf5_sequence
    path "${sample_id}_${ref_name}_ORF7_from_assembly.fasta", emit: orf7_sequence

    script:
    """
    set -e
    set -x

    echo "Extracting ORF sequences from assembly for sample ${sample_id} using ${ref_name} reference"

    # Define ORF coordinates based on reference type
    if [[ "${ref_name}" == "prrsV1" ]]; then
        ORF5_START=13494
        ORF5_END=14099
        ORF7_START=14598
        ORF7_END=14984
    elif [[ "${ref_name}" == "prrsV2" ]]; then
        ORF5_START=13788
        ORF5_END=14390
        ORF7_START=14889
        ORF7_END=15260
    else
        echo "Unknown reference type: ${ref_name}"
        exit 1
    fi

    echo "Using coordinates for ${ref_name}: ORF5=\$ORF5_START-\$ORF5_END, ORF7=\$ORF7_START-\$ORF7_END"

    # Extract reference ORF sequences for BLAST database
    samtools faidx ${reference}
    REF_NAME=\$(grep '^>' ${reference} | head -n 1 | cut -d ' ' -f 1 | sed 's/>//')
    
    # Extract ORF5 and ORF7 from reference
    samtools faidx ${reference} \${REF_NAME}:\$ORF5_START-\$ORF5_END > ref_orf5.fasta
    samtools faidx ${reference} \${REF_NAME}:\$ORF7_START-\$ORF7_END > ref_orf7.fasta

    # Create BLAST databases for ORF sequences
    makeblastdb -in ref_orf5.fasta -dbtype nucl -out orf5_db
    makeblastdb -in ref_orf7.fasta -dbtype nucl -out orf7_db

    # Function to find best matching contig for an ORF
    find_best_orf_match() {
        local orf_name=\$1
        local ref_orf=\$2
        local db_name=\$3
        
        echo "Searching for \$orf_name in assembly..."
        
        # BLAST assembly contigs against ORF reference
        blastn -query ${assembly} -db \$db_name -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" -max_target_seqs 1 -max_hsps 1 > \${orf_name}_blast_results.txt
        
        if [ -s "\${orf_name}_blast_results.txt" ]; then
            # Get the best hit
            best_hit=\$(head -n 1 \${orf_name}_blast_results.txt)
            contig_id=\$(echo "\$best_hit" | cut -f1)
            qstart=\$(echo "\$best_hit" | cut -f5)
            qend=\$(echo "\$best_hit" | cut -f6)
            pident=\$(echo "\$best_hit" | cut -f3)
            
            echo "Best \$orf_name match: contig \$contig_id, positions \$qstart-\$qend, identity: \$pident%"
            
            # Extract the ORF region from the contig
            if [ \$qstart -le \$qend ]; then
                samtools faidx ${assembly} \${contig_id}:\$qstart-\$qend > ${sample_id}_${ref_name}_\${orf_name}_from_assembly.fasta
            else
                # Reverse complement if needed
                samtools faidx ${assembly} \${contig_id}:\$qend-\$qstart > temp_\${orf_name}.fasta
                seqtk seq -r temp_\${orf_name}.fasta > ${sample_id}_${ref_name}_\${orf_name}_from_assembly.fasta
                rm temp_\${orf_name}.fasta
            fi
            
            # Rename sequence header
            sed -i "1s/.*/>${sample_id}_${ref_name}_\${orf_name}_from_assembly/" ${sample_id}_${ref_name}_\${orf_name}_from_assembly.fasta
        else
            echo "No significant match found for \$orf_name in assembly"
            touch ${sample_id}_${ref_name}_\${orf_name}_from_assembly.fasta
        fi
    }

    # Find ORF5 and ORF7 in assembly
    find_best_orf_match "ORF5" "ref_orf5.fasta" "orf5_db"
    find_best_orf_match "ORF7" "ref_orf7.fasta" "orf7_db"

    echo "ORF extraction completed for sample ${sample_id} with reference ${ref_name}"
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

process DOWNSAMPLE_BAM {
    conda "bioconda::samtools=1.15"
    
    input:
    tuple val(sample_id), val(reference_name), path(bam), path(bai)

    output:
    tuple val(sample_id), val(reference_name), path("${sample_id}_${reference_name}_downsampled.bam"), path("${sample_id}_${reference_name}_downsampled.bam.bai"), emit: downsampled_bam

    script:
    """
    set -e
    set -x

    # Count total reads in the BAM file
    total_reads=\$(samtools view -c ${bam})
    echo "Total reads in ${bam}: \$total_reads"

    # Set target reads to 100,000
    target_reads=100000

    if [ \$total_reads -le \$target_reads ]; then
        echo "BAM file has \$total_reads reads, which is <= \$target_reads. No downsampling needed."
        # Copy the original files with the correct output names
        cp ${bam} ${sample_id}_${reference_name}_downsampled.bam
        cp ${bai} ${sample_id}_${reference_name}_downsampled.bam.bai
    else
        echo "Downsampling from \$total_reads to \$target_reads reads"
        
        # Write a simple Python script to calculate the fraction
        cat > calc_fraction.py << 'EOF'
import sys
target = int(sys.argv[1])
total = int(sys.argv[2])
fraction = target / total
print(f"{fraction:.6f}")
EOF
        
        # Calculate the sampling fraction
        fraction=\$(python3 calc_fraction.py \$target_reads \$total_reads)
        echo "Sampling fraction: \$fraction"
        
        # Downsample the BAM file using samtools with seed for reproducibility
        samtools view -bs 42\$fraction ${bam} > ${sample_id}_${reference_name}_downsampled.bam
        
        # Index the downsampled BAM file
        samtools index ${sample_id}_${reference_name}_downsampled.bam
        
        # Report final read count
        final_reads=\$(samtools view -c ${sample_id}_${reference_name}_downsampled.bam)
        echo "Final read count after downsampling: \$final_reads"
    fi

    # Verify that both output files exist
    if [ ! -f "${sample_id}_${reference_name}_downsampled.bam" ]; then
        echo "ERROR: Output BAM file not created!"
        exit 1
    fi
    
    if [ ! -f "${sample_id}_${reference_name}_downsampled.bam.bai" ]; then
        echo "ERROR: Output BAM index file not created!"
        exit 1
    fi
    
    echo "Successfully created downsampled BAM and index files"
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
    tuple val(sample_id), val(ref_name), path("${sample_id}_${ref_name}_ORF5_consensus.fasta"), path("${sample_id}_${ref_name}_ORF7_consensus.fasta"), emit: consensus
    path "${sample_id}_${ref_name}_ORF5_variants.vcf.gz", emit: orf5_variants
    path "${sample_id}_${ref_name}_ORF7_variants.vcf.gz", emit: orf7_variants

    script:
    """
    set -e
    set -x

    # Get reference sequence name
    REF_NAME=\$(grep '^>' ${reference} | head -n 1 | cut -d ' ' -f 1 | sed 's/>//')
    echo "Reference sequence name: \$REF_NAME"

    # Define ORF coordinates based on reference type
    if [[ "${ref_name}" == "prrsV1" ]]; then
        ORF5_COORDS="13494-14099"
        ORF7_COORDS="14598-14984"
    elif [[ "${ref_name}" == "prrsV2" ]]; then
        ORF5_COORDS="13788-14390"
        ORF7_COORDS="14889-15260"
    else
        echo "Unknown reference type: ${ref_name}"
        exit 1
    fi

    echo "Using coordinates for ${ref_name}: ORF5=\$ORF5_COORDS, ORF7=\$ORF7_COORDS"

    # Function to create consensus for a specific region
    create_consensus_region() {
        local region=\$1
        local orf_name=\$2
        local coords=\$3
        
        echo "Processing \$orf_name region: \$region"
        
        # Check if there are any reads in the specified region
        read_count=\$(samtools view ${bam} \$region | wc -l)
        echo "Read count in \$region: \$read_count"

        if [ \$read_count -eq 0 ]; then
            echo "No reads found in \$region for sample ${sample_id} with reference ${ref_name}. Creating empty consensus."
            touch ${sample_id}_${ref_name}_\${orf_name}_consensus.fasta
            touch ${sample_id}_${ref_name}_\${orf_name}_variants.vcf.gz
        else
            # Generate VCF file for the specific region with minimum variant frequency of 10%
            bcftools mpileup -f ${reference} -q 20 -Q 20 -r \$region ${bam} | \
            bcftools call -mv --ploidy 1 -V indels -Oz -o ${sample_id}_${ref_name}_\${orf_name}_variants.vcf.gz

            # Index the VCF file
            bcftools index ${sample_id}_${ref_name}_\${orf_name}_variants.vcf.gz

            # Filter variants to keep only those with >= 10% frequency
            bcftools filter -i 'INFO/DP4[2]+INFO/DP4[3] >= 0.1*(INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])' \
                ${sample_id}_${ref_name}_\${orf_name}_variants.vcf.gz -Oz -o ${sample_id}_${ref_name}_\${orf_name}_filtered_variants.vcf.gz

            # Index the filtered VCF file
            bcftools index ${sample_id}_${ref_name}_\${orf_name}_filtered_variants.vcf.gz

            # Create consensus sequence with IUPAC codes for ambiguous positions (>= 10% frequency)
            bcftools consensus -f ${reference} --haplotype I -o ${sample_id}_${ref_name}_\${orf_name}_full_consensus.fasta ${sample_id}_${ref_name}_\${orf_name}_filtered_variants.vcf.gz

            # Extract the region of interest
            samtools faidx ${sample_id}_${ref_name}_\${orf_name}_full_consensus.fasta \$region > ${sample_id}_${ref_name}_\${orf_name}_consensus.fasta

            # Rename the sequence in the FASTA file to include sample and ORF info
            sed -i "1s/.*/>${sample_id}_${ref_name}_\${orf_name}/" ${sample_id}_${ref_name}_\${orf_name}_consensus.fasta

            # Clean up intermediate files
            rm ${sample_id}_${ref_name}_\${orf_name}_full_consensus.fasta*
        fi
    }

    # Create consensus for ORF5
    ORF5_REGION="\${REF_NAME}:\$ORF5_COORDS"
    create_consensus_region \$ORF5_REGION "ORF5" \$ORF5_COORDS

    # Create consensus for ORF7
    ORF7_REGION="\${REF_NAME}:\$ORF7_COORDS"
    create_consensus_region \$ORF7_REGION "ORF7" \$ORF7_COORDS

    echo "Consensus creation completed for sample ${sample_id} with reference ${ref_name}."
    echo "ORF5 consensus: ${sample_id}_${ref_name}_ORF5_consensus.fasta"
    echo "ORF7 consensus: ${sample_id}_${ref_name}_ORF7_consensus.fasta"
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
