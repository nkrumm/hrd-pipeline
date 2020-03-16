#!/usr/bin/env nextflow
def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --input sample.tsv -profile docker
""".stripIndent()
}
params.genome = 'GRCh37'
params.assay = 'GLTv3'

// Read in the genome fasta and index
ref_fasta = file(params.genomes[params.genome].ref_fasta, checkIfExists: true)
ref_fasta_fai = file(params.genomes[params.genome].ref_fasta_fai, checkIfExists: true)
gatk_mills = file(params.genomes[params.genome].gatk_mills, checkIfExists: true)
gatk_mills_index = file(params.genomes[params.genome].gatk_mills_index, checkIfExists: true)
gatk_1kg = file(params.genomes[params.genome].gatk_1kg, checkIfExists: true)
gatk_1kg_index = file(params.genomes[params.genome].gatk_1kg_index, checkIfExists: true)

bwa_index = Channel.fromPath(params.genomes[params.genome].bwa_index, checkIfExists: true).collect()
ref_fasta_dict = Channel.fromPath(params.genomes[params.genome].ref_fasta_dict, checkIfExists: true).collect()
//picard_intervals = file(params.assays[params.assay].picard_intervals, checkIfExists: true)
bed_file = file(params.assays[params.assay].bed_file, checkIfExists: true)

fastq_pair_ch = Channel.fromFilePairs(params.input + '/*{1,2}.fastq.gz', flat: true, checkIfExists: true)

process bwa_mem {
    // Align fastqs
    label 'bwa'

    tag "${sample_id}"

    input:
        path ref_fasta
        path bwa_index
        tuple val(sample_id), file(fastq1), file(fastq2) from fastq_pair_ch

    output:
        tuple val(sample_id), file("${sample_id}.sam") into align_ch

    cpus 16

    //publishDir params.output, overwrite: true

    // -Y use soft clipping for supplimentary alignment
    // -K process INT input bases in each batch regardless of nThreads (for reproducibility)
    // -C append FASTA/FASTQ comment to SAM output
    script:
    """ 
    bwa mem -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:NA\\tSM:${sample_id}\\t" -K 100000000 -t ${task.cpus} ${ref_fasta} ${fastq1} ${fastq2} > ${sample_id}.sam
    """
}

process samtools_view_sort {
    label 'samtools'

    tag "${sample_id}"

    input:
        tuple val(sample_id), file(sam_file) from align_ch

    output:
        tuple val(sample_id), file("${sample_id}.bam") into raw_bams

    cpus 16

    //publishDir params.output, overwrite: true

    script:
    """
    samtools view -hb $sam_file | samtools sort - -o ${sample_id}.bam
    """
}

process picard_remove_duplicates {
    label 'picard'

    tag "${sample_id}"

    input:
        tuple val(sample_id), file(bam_file) from raw_bams

    output:
        tuple val(sample_id), file("${sample_id}.bam") into rmdup_bams_recal, rmdup_bams

    cpus 16

    //publishDir params.output, overwrite: true

    script:
    """
    java -Xmx4g -jar /usr/picard/picard.jar MarkDuplicates INPUT=${bam_file} OUTPUT=${sample_id}.bam METRICS_FILE=${sample_id}.quality_metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true 2> picard_rmdupes.log
    """
}

process gatk_bqsr {
    label 'gatk'

    tag "${sample_id}"

    input:
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path gatk_mills
        path gatk_mills_index
        path gatk_1kg
        path gatk_1kg_index
        tuple val(sample_id), file(bam_file) from rmdup_bams_recal

    output:
        tuple val(sample_id), file("${sample_id}.recal_table") into bqsr_recal_tables

    cpus 16

    //publishDir params.output, overwrite: true

    script:
    """
    gatk --java-options "-Xmx4G" BaseRecalibrator --reference ${ref_fasta} --input ${bam_file} --known-sites ${gatk_mills} --known-sites ${gatk_1kg} --output ${sample_id}.recal_table 2> gatk_bqsr.log
    """
}

process gatk_apply_bqsr {
    label 'gatk'

    tag "${sample_id}"

    bqsr_apply = rmdup_bams.join(bqsr_recal_tables)

    input:
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        tuple val(sample_id), file(bam_file), file(recal_table) from bqsr_apply

    output:
        tuple val(sample_id), file("${sample_id}.bqsr.bam") into bqsr_bams

    cpus 16

    //publishDir params.output, overwrite: true

    script:
    """
    gatk ApplyBQSR --reference ${ref_fasta} --input ${bam_file} --bqsr-recal-file ${recal_table} --output ${sample_id}.bqsr.bam 2> gatk_apply_bqsr.log
    """
}

process samtools_final_bam {
    label 'samtools'

    tag "${sample_id}"

    input:
        tuple val(sample_id), file(bqsr_bam) from bqsr_bams

    output:
        tuple val(sample_id), file("${sample_id}.final.bam") into final_bams

    cpus 16

    publishDir params.output, overwrite: true

    script:
    """
    samtools sort ${bqsr_bam} -o ${sample_id}.final.bam && samtools index ${sample_id}.final.bam
    """
}

process samtools_mpileup {
    label 'samtools'

    tag "${sample_id}"

    input:
        path ref_fasta
        path bed_file
        tuple val(sample_id), file(final_bam) from final_bams

    output:
        tuple val(sample_id), file("${sample_id}.mpileup") into mpileup

    cpus 16

    publishDir params.output, overwrite: true

    shell:
    '''
    samtools mpileup -f !{ref_fasta} -d 1000000 -A -E -l !{bed_file} !{final_bam} 2> mpileup.log | awk '{if($4 > 0) print $0}' > !{sample_id}.mpileup || rm !{sample_id}.mpileup
    '''
}
