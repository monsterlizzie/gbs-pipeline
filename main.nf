nextflow.enable.dsl=2
// import modules for process
include { FILE_VALIDATION; PREPROCESS; READ_QC } from "$projectDir/modules/preprocess"
include { ASSEMBLY_UNICYCLER; ASSEMBLY_SHOVILL; ASSEMBLY_ASSESS; ASSEMBLY_QC; ASSEMBLY_QC_FALLBACK } from "$projectDir/modules/assembly"
include { GET_REF_GENOME_BWA_DB; MAPPING; SAM_TO_SORTED_BAM; SNP_CALL; HET_SNP_COUNT; MAPPING_QC; MAPPING_QC_FALLBACK} from "$projectDir/modules/mapping"
include { GET_KRAKEN2_DB; TAXONOMY; BRACKEN;TAXONOMY_QC; TAXONOMY_QC_FALLBACK } from "$projectDir/modules/taxonomy"
include { OVERALL_QC } from "$projectDir/modules/overall_qc"
include { GENERATE_SAMPLE_REPORT; GENERATE_OVERALL_REPORT } from "$projectDir/modules/output"
include { serotyping } from './modules/serotyping.nf'
include { srst2_for_res_typing; split_target_RES_seq_from_sam_file; split_target_RES_sequences; freebayes } from './modules/res_alignments.nf'
include { res_typer } from './modules/res_typer.nf'
include { surface_typer } from './modules/surface_typer.nf'
include { getmlst_for_srst2; srst2_for_mlst; get_mlst_allele_and_pileup} from './modules/mlst.nf'
include { get_pbp_genes; get_pbp_alleles } from './modules/pbp_typer.nf'
include { finalise_sero_res_results; finalise_surface_typer_results; finalise_pbp_existing_allele_results; combine_results } from './modules/combine.nf'
include { get_version } from './modules/version.nf'

// Add this utility process to ensure 'databases/' exists
process INIT_DB_DIR {
    label 'bash_container'
    label 'farm_local'
    publishDir "${params.db}"
   

    output:
    val "${params.db}", emit: db_dir
    path "do_not_modify", emit:dummy


    script:
    """
    mkdir do_not_modify
    """
}

// Main pipeline workflow
workflow {

    main:

    INIT_DB_DIR()
    db_dir_ch = INIT_DB_DIR.out.db_dir
    
    // Get path and prefix of Reference Genome BWA Database, generate from assembly if necessary
    GET_REF_GENOME_BWA_DB(params.ref_genome, db_dir_ch)

    // Get path to Kraken2 Database, download if necessary
    GET_KRAKEN2_DB(params.kraken2_db_remote, db_dir_ch)

    // Bracken too, GET_KRAKEN2_DB(params.kraken2_db_remote, params.db)
 
// Get read pairs into Channel raw_read_pairs_ch
    raw_read_pairs_ch = Channel.fromFilePairs("$params.reads/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}", checkIfExists: true)

    // Basic input files validation
    // Output into Channel FILE_VALIDATION.out.result
    FILE_VALIDATION(raw_read_pairs_ch)

    // From Channel raw_read_pairs_ch, only output valid reads of samples based on Channel FILE_VALIDATION.out.result
    VALID_READS_ch = FILE_VALIDATION.out.result.join(raw_read_pairs_ch, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    // Preprocess valid read pairs
    // Output into Channels PREPROCESS.out.processed_reads & PREPROCESS.out.json
    PREPROCESS(VALID_READS_ch)

    // From Channel PREPROCESS.out.json, provide Read QC status
    // Output into Channels READ_QC.out.bases, READ_QC.out.result, READ_QC.out.report
    READ_QC(PREPROCESS.out.json, params.length_low, params.depth)

    // From Channel PREPROCESS.out.processed_reads, only output reads of samples passed Read QC based on Channel READ_QC.out.result
    READ_QC_PASSED_READS_ch = READ_QC.out.result.join(PREPROCESS.out.processed_reads, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    
    // From Channel READ_QC_PASSED_READS_ch, assemble the preprocess read pairs
    // Output into Channel ASSEMBLY_ch, and hardlink (default) the assemblies to $params.output directory
    switch (params.assembler) {
        case 'shovill':
            ASSEMBLY_ch = ASSEMBLY_SHOVILL(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break

        case 'unicycler':
            ASSEMBLY_ch = ASSEMBLY_UNICYCLER(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break
    }

    // From Channel ASSEMBLY_ch, assess assembly quality
    // Output into Channel ASSEMBLY_ASSESS.out.report
    ASSEMBLY_ASSESS(ASSEMBLY_ch)

    // From Channel ASSEMBLY_ASSESS.out.report and Channel READ_QC.out.bases, provide Assembly QC status
    // Output into Channels ASSEMBLY_QC.out.result & ASSEMBLY_QC.out.report
    ASSEMBLY_QC(
        ASSEMBLY_ASSESS.out.report
        .join(READ_QC.out.bases, failOnDuplicate: true),
        params.contigs,
        params.length_low,
        params.length_high,
        params.depth
    )

    // Provide fallback Assembly QC for samples that failed READ_QC
    READ_QC_FAIL_SAMPLE_ID_ch = READ_QC.out.result.filter { it[1] == "FAIL" }.map { it[0] }
    ASSEMBLY_QC_FALLBACK(READ_QC_FAIL_SAMPLE_ID_ch)

    assembly_qc_all = ASSEMBLY_QC.out.result.mix(ASSEMBLY_QC_FALLBACK.out.result)
    assembly_qc_report_all = ASSEMBLY_QC.out.report.mix(ASSEMBLY_QC_FALLBACK.out.report)

    // From Channel READ_QC_PASSED_READS_ch map reads to reference
    // Output into Channel MAPPING.out.sam
    MAPPING(GET_REF_GENOME_BWA_DB.out.path, GET_REF_GENOME_BWA_DB.out.prefix, READ_QC_PASSED_READS_ch)

    // From Channel MAPPING.out.sam, Convert SAM into sorted BAM and calculate reference coverage
    // Output into Channels SAM_TO_SORTED_BAM.out.sorted_bam and SAM_TO_SORTED_BAM.out.ref_coverage
    SAM_TO_SORTED_BAM(MAPPING.out.sam, params.lite)

    // From Channel SAM_TO_SORTED_BAM.out.sorted_bam calculates non-cluster Het-SNP site count
    // Output into Channel HET_SNP_COUNT.out.result
    SNP_CALL(params.ref_genome, SAM_TO_SORTED_BAM.out.sorted_bam, params.lite)
    HET_SNP_COUNT(SNP_CALL.out.vcf)

    // Merge Channels SAM_TO_SORTED_BAM.out.ref_coverage & HET_SNP_COUNT.out.result to provide Mapping QC Status
    // Output into Channels MAPPING_QC.out.result & MAPPING_QC.out.report
    MAPPING_QC(
    SAM_TO_SORTED_BAM.out.ref_coverage
    .join(HET_SNP_COUNT.out.result, failOnDuplicate: true),
    params.ref_coverage,
    params.het_snp_site
)

    // Provide fallback Mapping QC for samples that failed READ_QC
    READ_QC_FAIL_ch = READ_QC.out.result.filter { it[1] == "FAIL" }
    MAPPING_QC_FALLBACK(READ_QC_FAIL_ch, params.ref_coverage, params.het_snp_site)

    // Now safely mix both outputs
    mapping_qc_all = MAPPING_QC.out.result.mix(MAPPING_QC_FALLBACK.out.result)
    mapping_qc_report_all = MAPPING_QC.out.report.mix(MAPPING_QC_FALLBACK.out.report)


    // From Channel READ_QC_PASSED_READS_ch assess Streptococcus agalactiae percentage in reads
    // Output into Channel TAXONOMY.out.report
    TAXONOMY(GET_KRAKEN2_DB.out.path, params.kraken2_memory_mapping, READ_QC_PASSED_READS_ch)
    //Bracken too, estimate the abundance of Streptococcus agalactiae
    BRACKEN(
        GET_KRAKEN2_DB.out.path,
        TAXONOMY.out.report,
        params.read_len,
        params.classification_level,
        params.threshold
    )
    // From Channel TAXONOMY.out.report, provide taxonomy QC status
    // Output into Channels TAXONOMY_QC.out.result & TAXONOMY_QC.out.report
    TAXONOMY_QC(BRACKEN.out.bracken_report, params.sagalactiae_percentage, params.top_non_agalactiae_species_percentage)

    TAXONOMY_QC_FALLBACK(READ_QC_FAIL_SAMPLE_ID_ch)
    taxonomy_qc_all = TAXONOMY_QC.out.result.mix(TAXONOMY_QC_FALLBACK.out.result)
    taxonomy_qc_report_all = TAXONOMY_QC.out.report.mix(TAXONOMY_QC_FALLBACK.out.report)


    // Merge Channels FILE_VALIDATION.out.result & READ_QC.out.result & ASSEMBLY_QC.out.result & MAPPING_QC.out.result & TAXONOMY_QC.out.result to provide Overall QC Status
    // Output into Channel OVERALL_QC.out.result & OVERALL_QC.out.report
    OVERALL_QC(
        raw_read_pairs_ch.map{ it[0] }
        .join(FILE_VALIDATION.out.result, failOnDuplicate: true, remainder: true)
        .join(READ_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(assembly_qc_all, failOnDuplicate: true, remainder: true)
        .join(mapping_qc_all, failOnDuplicate: true, remainder: true)
        .join(taxonomy_qc_all, failOnDuplicate: true, remainder: true)
    )

    // From Channel READ_QC_PASSED_READS_ch, only output reads of samples passed overall QC based on Channel OVERALL_QC.out.result
    OVERALL_QC_PASSED_READS_ch = OVERALL_QC.out.result.join(READ_QC_PASSED_READS_ch, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    // Extract only paired reads (R1 and R2) from OVERALL_QC_PASSED_READS_ch
    OVERALL_QC_PASSED_PAIRED_READS_ch = OVERALL_QC_PASSED_READS_ch.map { id, r1, r2, unpaired -> tuple(id, r1, r2) }
    
    // From Channel ASSEMBLY_ch, only output assemblies of samples passed overall QC based on Channel OVERALL_QC.out.result
    OVERALL_QC_PASSED_ASSEMBLIES_ch = OVERALL_QC.out.result.join(ASSEMBLY_ch, failOnDuplicate: true)
                            .filter { it[1] == 'PASS' }
                            .map { it[0, 2..-1] }


    // Generate sample reports by merging outputs from all result-generating modules
    GENERATE_SAMPLE_REPORT(
        raw_read_pairs_ch.map{ it[0] }
        .join(READ_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(assembly_qc_report_all, failOnDuplicate: true, remainder: true)
        .join(mapping_qc_report_all, failOnDuplicate: true, remainder: true)
        .join(taxonomy_qc_report_all, failOnDuplicate: true, remainder: true)
        .join(OVERALL_QC.out.report, failOnDuplicate: true, failOnMismatch: true)
        .map { row ->
            def sample_id = row[0]
            def report_paths = row[1..-1].findAll { it != null && it != "NA" }
            [sample_id, report_paths]
        }
    )



        // Check if reads specified
    if (params.run_sero_res | params.run_mlst | params.run_surfacetyper){
        if (params.reads == ""){
            println("Please specify reads with --reads.")
            println("Print help with nextflow main.nf --help")
            System.exit(1)
        }
        // Create read pairs channel
        Channel.fromFilePairs( params.reads, checkIfExists: true )
            .set { read_pairs_ch }
    }

    // Check if params.output specified
    if (params.output == ""){
        println("Please specify the results directory with --params.output.")
        println("Print help with nextflow main.nf --help")
        System.exit(1)
    }

    if (!params.run_sero_res && !params.run_surfacetyper && !params.run_mlst && !params.run_pbptyper){
        println("Please specify one or more pipelines to run.")
        println("Print help with nextflow main.nf --help")
        System.exit(1)
    }

    // Check parameters are within range
    if (params.gbs_res_min_coverage < 0 | params.gbs_res_min_coverage > 100){
        println("--gbs_res_min_coverage value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }

    if (params.gbs_res_max_divergence < 0 | params.gbs_res_max_divergence > 100){
        println("--gbs_res_max_divergence value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }

    other_res_min_coverage_list = params.other_res_min_coverage.toString().tokenize(' ')
    for (other_res_min_coverage in other_res_min_coverage_list){
        if (other_res_min_coverage.toDouble() < 0 | other_res_min_coverage.toDouble() > 100){
            println("--other_res_min_coverage value(s) not in range. Please specify a value between 0 and 100.")
            System.exit(1)
        }
    }

    other_res_max_divergence_list = params.other_res_max_divergence.toString().tokenize(' ')
    for (other_res_max_divergence in other_res_max_divergence_list){
        if (other_res_max_divergence.toDouble() < 0 | other_res_max_divergence.toDouble() > 100){
            println("--other_res_max_divergence value(s) not in range. Please specify a value between 0 and 100.")
            System.exit(1)
        }
    }

    if (params.restyper_min_read_depth < 0){
        println("--restyper_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }

    if (params.serotyper_min_read_depth < 0){
        println("--serotyper_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }

    if (params.mlst_min_coverage < 0 | params.mlst_min_coverage > 100){
        println("--mlst_min_coverage value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }

    if (params.mlst_min_read_depth < 0){
        println("--mlst_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }

    if (params.surfacetyper_min_coverage < 0 | params.surfacetyper_min_coverage > 100){
        println("--surfacetyper_min_coverage value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }

    if (params.surfacetyper_max_divergence < 0 | params.surfacetyper_max_divergence > 100){
        println("--surfacetyper_max_divergence value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }

    if (params.surfacetyper_min_read_depth < 0){
        println("--surfacetyper_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }

    
    // Output files
    params.sero_res_incidence_out         = "${params.output}/typer/serotype_res_incidence.txt"
    params.variants_out                   = "${params.output}/typer/gbs_res_variants.txt"
    params.alleles_variants_out          = "${params.output}/typer/drug_cat_alleles_variants.txt"
    params.res_accessions_out            = "${params.output}/typer/resfinder_accessions.txt"
    params.existing_pbp_alleles_out      = "${params.output}/typer/existing_pbp_alleles.txt"
    params.surface_protein_incidence_out = "${params.output}/typer/surface_protein_incidence.txt"
    params.surface_protein_variants_out  = "${params.output}/typer/surface_protein_variants.txt"
    params.existing_mlst_alleles_out     = "${params.output}/typer/existing_sequence_types.txt"
    params.new_mlst_alleles_status       = "${params.output}/typer/new_mlst_alleles.log"
    params.gbs_typer_report              = "${params.output}/typer/gbs_typer_report.txt"


   


    workflow GBS_RES {

    take:
        reads

    main:

        gbs_res_typer_db = file(params.gbs_res_typer_db, checkIfExists: true)
        gbs_res_targets_db = file(params.gbs_res_targets_db, checkIfExists: true)

        // Split GBS target sequences from GBS resistance database into separate FASTA files per sequence
        split_target_RES_sequences(gbs_res_typer_db, gbs_res_targets_db)

        // Map genomes to GBS resistance database using SRST2
        srst2_for_res_typing(reads, gbs_res_typer_db, params.gbs_res_min_coverage, params.gbs_res_max_divergence)
        fullgenes = srst2_for_res_typing.out.fullgenes

        // Split sam file for each GBS target sequence
        split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, gbs_res_targets_db)

        // Get consensus sequence using freebayes
        freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out)
        consensus = freebayes.out.consensus

    emit:
        fullgenes
        consensus
}

// Resistance mapping with the other resistance databases
workflow OTHER_RES {

    take:
        reads

    main:
        other_res_db = file(params.other_res_db, checkIfExists: true)
        // Map genomes to resistance database using SRST2
        srst2_for_res_typing(reads, other_res_db, params.other_res_min_coverage, params.other_res_max_divergence)
        fullgenes = srst2_for_res_typing.out.fullgenes

    emit:
        fullgenes
}

// MLST pipeline
workflow MLST {

    take:
        reads

    main:
        // Get MLST database for all downstream processes
        getmlst_for_srst2()

        // Run SRST2 MLST
        srst2_for_mlst(getmlst_for_srst2.out.getmlst_results, reads, params.mlst_min_coverage)

        // Get new consensus allele and pileup data
        get_mlst_allele_and_pileup(srst2_for_mlst.out.bam_and_srst2_results, params.mlst_min_read_depth)

        // Collect outputs
        new_alleles = get_mlst_allele_and_pileup.out.new_alleles
        pileup = get_mlst_allele_and_pileup.out.pileup
        existing_alleles = get_mlst_allele_and_pileup.out.existing_alleles
        status = get_mlst_allele_and_pileup.out.new_alleles_status
        srst2_results = srst2_for_mlst.out.srst2_results

    emit:
        new_alleles
        pileup
        existing_alleles
        status
        srst2_results
}

// PBP-1A allele typing pipeline
workflow PBP1A {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS1A-1', file(params.gbs_blactam_1A_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${params.output}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// PBP-2B allele typing pipeline
workflow PBP2B {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS2B-1', file(params.gbs_blactam_2B_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${params.output}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// PBP-2X allele typing pipeline
workflow PBP2X {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS2X-1', file(params.gbs_blactam_2X_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${params.output}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

 if (params.run_sero_res){

        // Serotyping Process
        serotyping(OVERALL_QC_PASSED_PAIRED_READS_ch, file(params.sero_gene_db, checkIfExists: true), params.serotyper_min_read_depth)

        // Resistance Mapping Workflows
        GBS_RES(OVERALL_QC_PASSED_PAIRED_READS_ch)
        OTHER_RES(OVERALL_QC_PASSED_PAIRED_READS_ch)

        // Once GBS or both resistance workflows are complete, trigger resistance typing
        GBS_RES.out.fullgenes
        .join(GBS_RES.out.consensus)
        .join(OTHER_RES.out.fullgenes)
        .set { res_files_ch }

        res_typer(res_files_ch, params.restyper_min_read_depth, file(params.config, checkIfExists: true))

        // Combine serotype and resistance type results for each sample
        sero_res_ch = serotyping.out.join(res_typer.out.res_out)

        finalise_sero_res_results(sero_res_ch, file(params.config, checkIfExists: true))

        // Combine samples and output results files
        finalise_sero_res_results.out.sero_res_incidence
            .collectFile(name: file("${params.output}/${params.sero_res_incidence_out}"), keepHeader: true)

        finalise_sero_res_results.out.res_alleles_variants
            .collectFile(name: file("${params.output}/${params.alleles_variants_out}"), keepHeader: true)

        finalise_sero_res_results.out.res_variants
            .collectFile(name: file("${params.output}/${params.variants_out}"), keepHeader: true)

        res_typer.out.res_accessions
            .collectFile(name: file("${params.output}/${params.res_accessions_out}"))
    }

    // MLST
    if (params.run_mlst){

        MLST(OVERALL_QC_PASSED_PAIRED_READS_ch)
        MLST.out.new_alleles.subscribe { it ->
            it.copyTo(file("${params.output}"))
        }
        MLST.out.pileup.subscribe { it ->
            it.copyTo(file("${params.output}"))
        }
        MLST.out.existing_alleles
            .collectFile(name: file("${params.output}/${params.existing_mlst_alleles_out}"), keepHeader: true, sort: true)
        MLST.out.status
            .collectFile(name: file("${params.output}/${params.new_mlst_alleles_status}"), keepHeader: false, sort: true)

    }

    // Surface Typing Process
    if (params.run_surfacetyper){

        surface_typer(OVERALL_QC_PASSED_PAIRED_READS_ch, file(params.gbs_surface_typer_db, checkIfExists: true),
            params.surfacetyper_min_read_depth, params.surfacetyper_min_coverage,
            params.surfacetyper_max_divergence)

        finalise_surface_typer_results(surface_typer.out, file(params.config, checkIfExists: true))

        // Combine results for surface typing
        finalise_surface_typer_results.out.surface_protein_incidence
            .collectFile(name: file("${params.output}/${params.surface_protein_incidence_out}"), keepHeader: true)
        finalise_surface_typer_results.out.surface_protein_variants
            .collectFile(name: file("${params.output}/${params.surface_protein_variants_out}"), keepHeader: true)

    }

    // PBP Typer
    if (params.run_pbptyper){

        // Check if contigs specified
        if (params.contigs == ""){
            println("Please specify contigs with --contigs.")
            println("Print help with --contigs")
            System.exit(1)
        }

        contig_paths = OVERALL_QC_PASSED_ASSEMBLIES_ch
            .fromPath(params.contigs, checkIfExists: true)
            .map { file -> tuple(file.baseName, file) }

        get_pbp_genes(contig_paths, file(params.gbs_blactam_db, checkIfExists: true), params.pbp_frac_align_threshold, params.pbp_frac_identity_threshold)

        // Get PBP existing and new alleles
        PBP1A(get_pbp_genes.out)
        PBP2B(get_pbp_genes.out)
        PBP2X(get_pbp_genes.out)

        PBP1A.out
        .concat(PBP2B.out, PBP2X.out)
        .set { PBP_all }

        PBP_all
            .collectFile(name: file("${params.output}/${params.existing_pbp_alleles_out}"), keepHeader: true, sort: true)
    }

    // Combine serotype, resistance, allelic profile, surface typer and GBS resistance variants
    if (params.run_sero_res & params.run_surfacetyper & params.run_mlst){

        // Get version of pipeline
        get_version()
        version_ch = get_version.out

        // Combine serotype and resistance type results for each sample
        combined_ch = serotyping.out
            .join(res_typer.out.res_out)
            .join(surface_typer.out)
            .join(MLST.out.srst2_results)

        combine_results(combined_ch, file(params.config, checkIfExists: true), version_ch)

        combine_results.out
            .collectFile(name: file("${params.output}/${params.gbs_typer_report}"), keepHeader: true, sort: true)
    }

    //GENERATE_OVERALL_REPORT(GENERATE_SAMPLE_REPORT.out.report.collect()) for later combining qc and tyoer reports 

}