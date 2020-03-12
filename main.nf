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
params.assay = 'MONCv1'

// Read in the genome fasta and index
ref_fasta = file(params.genomes[params.genome].ref_fasta, checkIfExists: true)
ref_fasta_fai = file(params.genomes[params.genome].ref_fasta_fai, checkIfExists: true)

bwa_index = Channel.fromPath(params.genomes[params.genome].bwa_index, checkIfExists: true).collect()
ref_fasta_dict = Channel.fromPath(params.genomes[params.genome].ref_fasta_dict, checkIfExists: true).collect()
//picard_intervals = file(params.assays[params.assay].picard_intervals, checkIfExists: true)
bed_file = file(params.assays[params.assay].bed_file, checkIfExists: true)

fastq_pair_ch = Channel.fromFilePairs('test_data/*{1,2}.fastq.gz', flat: true, checkIfExists: true)

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

    publishDir params.output, overwrite: true

    // -Y use soft clipping for supplimentary alignment
    // -K process INT input bases in each batch regardless of nThreads (for reproducibility)
    // -C append FASTA/FASTQ comment to SAM output
    script:
    """ 
    bwa mem -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" -K 10000000 -t ${task.cpus} ${ref_fasta} ${fastq1} ${fastq2} > ${sample_id}.sam
    """
}

process samtools_view_sort {
    label 'samtools'

    tag "${sample_id}"

    input:
        tuple val(sample_id), file(sam_file) from align_ch

    output:
        tuple val(sample_id), file("${sample_id}.bam") into raw_bam

    cpus 16

    publishDir params.output, overwrite: true

    script:
    """
    samtools view -h -u $sam_file | samtools sort - -o ${sample_id}.bam
    """
}

process picard_remove_duplicates {
    label 'picard'

    tag "${sample_id}"

    input:
        tuple val(sample_id), file(bam_file) from raw_bam

    output:
        tuple val(sample_id), file("${sample_id}.bam") into rmdup_bam

    cpus 16

    publishDir params.output, overwrite: true

    script:
    """
    java -Xmx4g -jar /usr/picard/picard.jar MarkDuplicates INPUT=${bam_file} OUTPUT=${sample_id}.bam METRICS_FILE=${sample_id}.quality_metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true 2> picard_rmdupes.log
    """
}
