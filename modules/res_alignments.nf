process srst2_for_res_typing {
    label 'srst2'
    label 'farm_mid'
    
    input:
    tuple val(pair_id), path(reads) // ID and paired read files
    path db // File of resistance database file(s)
    val(min_coverage) // String of minimum coverage parameter(s) for SRST2
    val(max_divergence) // String of maximum coverage parameter(s) for SRST2

    output:
    tuple val(pair_id), path("${pair_id}*.bam"), emit: bam_files
    tuple val(pair_id), path("${pair_id}__fullgenes__${db_name}__results.txt"), emit: fullgenes

    script:
    db_name=db.getSimpleName()
    """

    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ${pair_id} --log --save_scores --min_coverage ${min_coverage} --max_divergence ${max_divergence} --gene_db ${db}

    touch ${pair_id}__fullgenes__${db_name}__results.txt
    """
}

process split_target_RES_sequences {

    input:
    path(fasta_file) // FASTA file of GBS target sequences
    path(targets_file) // Text file of GBS targets of interest

    output:
    path("CHECK_*")

    """
    get_targets_from_db.py -f ${fasta_file} -t ${targets_file} -o CHECK_

    # Clean
    unlink ${fasta_file}
    unlink ${targets_file}
    """
}

process split_target_RES_seq_from_sam_file {
    label 'farm_mid'
    
    input:
    tuple val(pair_id), path(bam_file) // ID and corresponding BAM file from mapping
    path(targets_file) // Text file of GBS targets of interest

    output:
    val(pair_id)
    path("*_*_${pair_id}*.bam")
    path("*_*_${pair_id}*.bai")

    """
    set +e
    samtools view -h ${bam_file} > \$(basename ${bam_file} .bam).sam
    get_targets_from_samfile.py -s \$(basename ${bam_file} .bam).sam -t ${targets_file} -i ${pair_id} -o CHECK_
    for check_sam_file in CHECK_*_${pair_id}*.sam; do
        samtools view -bS \${check_sam_file} > \$(basename \${check_sam_file} .sam).bam
        samtools index \$(basename \${check_sam_file} .sam).bam \$(basename \${check_sam_file} .sam).bai
    done

    touch dummy_dummy_${pair_id}_dummy.bam
    touch dummy_dummy_${pair_id}_dummy.bai

    # Clean directory
    unlink ${bam_file}
    unlink ${targets_file}
    """
}

process freebayes {
    label 'farm_mid'

    input:
    val(pair_id) // ID
    path(target_bam) // BAM file from a mapped target sequence of interest
    path(target_bai) // Corresponding BAM index file
    path(target_ref) // FASTA file of target sequence

    output:
    tuple val(pair_id), path("${pair_id}_consensus_seq.fna"), emit: consensus

    """
    for check_bam_file in CHECK_*_${pair_id}*.bam; do
        target=\$(echo \${check_bam_file} | sed 's/CHECK_//g' | sed 's/_${pair_id}.*//g')
        freebayes -q 20 -p 1 -f CHECK_\${target}_ref.fna \${check_bam_file} -v CHECK_\${target}_${pair_id}_seq.vcf
        bgzip CHECK_\${target}_${pair_id}_seq.vcf
        tabix -p vcf CHECK_\${target}_${pair_id}_seq.vcf.gz
        cat CHECK_\${target}_ref.fna | vcf-consensus CHECK_\${target}_${pair_id}_seq.vcf.gz >> ${pair_id}_consensus_seq.fna
        rm CHECK_\${target}_${pair_id}_seq.vcf.gz
        rm CHECK_\${target}_${pair_id}_seq.vcf.gz.tbi
        rm CHECK_\${target}_ref.fna.fai
    done

    touch ${pair_id}_consensus_seq.fna

    # Clean directory
    for check_bam_file in CHECK_*_${pair_id}*.bam; do
        target=\$(echo \${check_bam_file} | sed 's/CHECK_//g' | sed 's/_${pair_id}.*//g')
        unlink \${check_bam_file}
        unlink CHECK_\${target}_${pair_id}_seq.bai
        unlink CHECK_\${target}_ref.fna
    done

    mkdir output
    mv ${pair_id}_consensus_seq.fna output
    find . -maxdepth 1 -type f -delete
    mv output/${pair_id}_consensus_seq.fna .
    rm -d output
    """
}
