#!/usr/bin/env nextflow
def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    ./nextflow run main.nf --tumor samples/ --normal samples/ --gc hg19.gc.wig --cen hg19.centromere.txt -profile docker
""".stripIndent()
}
params.genome = 'GRCh37'

// Read in params.genome files
ref_fasta = file(params.genomes[params.genome].ref_fasta, checkIfExists: true)
ref_index = Channel.fromPath(params.genomes[params.genome].ref_index, checkIfExists: true).collect()
gatk_mills = file(params.genomes[params.genome].gatk_mills, checkIfExists: true)
gatk_mills_index = file(params.genomes[params.genome].gatk_mills_index, checkIfExists: true)
gatk_1kg = file(params.genomes[params.genome].gatk_1kg, checkIfExists: true)
gatk_1kg_index = file(params.genomes[params.genome].gatk_1kg_index, checkIfExists: true)

// Read in command line input files
gc_window = file(params.gc_window, checkIfExists: true)
centromere_file = file(params.centromere, checkIfExists: true)
Channel.fromFilePairs(params.normal + '/*{1,2}.fastq.gz', flat: true, checkIfExists: true)
                       .map { tuple( it[0], "normal", it[1], it[2] ) }
                       .set {normal_samples}

Channel.fromFilePairs(params.tumor + '/*{1,2}.fastq.gz', flat: true, checkIfExists: true)
                       .map { tuple( it[0], "tumor", it[1], it[2] ) }
                       .set {tumor_samples}

fastqs = normal_samples.mix(tumor_samples)


process alignment {
    // Align fastqs, sort and index
    label 'alignment'

    tag "${sample_id}-${sample_type}"

    input:
        path ref_fasta
        path ref_index
        tuple val(sample_id), val(sample_type), file(fastq1), file(fastq2) from fastqs

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.raw.bam") into raw_bams

    // publishDir params.output, mode: 'copy', overwrite: true

    // -K process INT input bases in each batch regardless of nThreads (for reproducibility)
    script:
    """ 
    bwa mem \
       -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:NA\\tSM:${sample_id}\\t" \
       -K 100000000 \
       -t ${task.cpus}  \
       ${ref_fasta} ${fastq1} ${fastq2} 2> log.txt \
     | samtools sort -t${task.cpus} -m4G - -o ${sample_id}.${sample_type}.raw.bam
     """
}


process picard_remove_duplicates {
    label 'picard'

    tag "${sample_id}-${sample_type}"

    input:
        tuple val(sample_id), val(sample_type), file(bam_file) from raw_bams

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.rmdup.bam"), file("${sample_id}.${sample_type}.rmdup.bai") into rmdup_bams

    //publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ -Dpicard.useLegacyParser=false \
    MarkDuplicates \
    -INPUT ${bam_file} \
    -OUTPUT ${sample_id}.${sample_type}.rmdup.bam \
    -METRICS_FILE ${sample_id}.${sample_type}.quality_metrics \
    -REMOVE_DUPLICATES true \
    -ASSUME_SORTED true \
    -VALIDATION_STRINGENCY SILENT \
    -CREATE_INDEX true 2> picard_rmdupes.log
    """
}

process gatk_bqsr {
    label 'gatk'

    tag "${sample_id}-${sample_type}"

    input:
        path ref_fasta
        path ref_index
        path gatk_mills
        path gatk_mills_index
        path gatk_1kg
        path gatk_1kg_index
        tuple val(sample_id), val(sample_type), file(bam_file), file(bam_bai) from rmdup_bams

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.bqsr.bam") into bqsr_bams

    //publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    BaseRecalibrator \
    --reference ${ref_fasta} \
    --input ${bam_file} \
    --known-sites ${gatk_mills} \
    --known-sites ${gatk_1kg} \
    --output ${sample_id}.${sample_type}.recal_table

    gatk ApplyBQSR \
    --reference ${ref_fasta} \
    --input ${bam_file} \
    --bqsr-recal-file ${sample_id}.${sample_type}.recal_table \
    --output ${sample_id}.${sample_type}.bqsr.bam
    """
}


process samtools_final_bam {
    label 'samtools'

    tag "${sample_id}-${sample_type}"

    input:
        tuple val(sample_id), val(sample_type), file(bqsr_bam) from bqsr_bams

    output:
        tuple val(sample_id), val(sample_type), file("${sample_id}.${sample_type}.final.bam") into final_bams
        file("*.bai")

    publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    samtools sort ${bqsr_bam} -o ${sample_id}.${sample_type}.final.bam

    samtools index ${sample_id}.${sample_type}.final.bam
    """
}

process samtools_mpileup {
    label 'samtools'

    tag "${sample_id}-${sample_type}"

    input:
        path ref_fasta
        tuple val(sample_id), val(sample_type), file(final_bam) from final_bams

    output:
        tuple val(sample_id), file("${sample_id}.${sample_type}.mpileup") into mpileups

    publishDir params.output, mode: 'copy', overwrite: true

    // -l ${bed_file}
    // -B --no-BAQ
    // -E --redo-BAQ
    script:
    """
    samtools mpileup -f ${ref_fasta} -d 1000000 -A -B ${final_bam} > ${sample_id}.${sample_type}.mpileup
    """
}

process sequenza_pileup2seqz {
    label 'sequenza'

    tag "${sample_id}"

    echo true

    paired_mpileups = mpileups
        .groupTuple(size: 2)
        .map{
            sample_id, mpileup_files -> tuple( sample_id, mpileup_files.sort{ it.getName() } )
            }

    input:
        path gc_window
        tuple val(sample_id), file(mpileup_files) from paired_mpileups
        // [sample_id, [normal.mpileup, tumor.mpileup] ]

    output:
        tuple val(sample_id), file("${sample_id}.seqz") into sequenza_seqz

    publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    sequenza-utils bam2seqz -gc ${gc_window} -p -n ${mpileup_files[0]} -t ${mpileup_files[1]} -o ${sample_id}.seqz
    """
}

process sequenza_seqz_binning {
    label 'sequenza'

    tag "${sample_id}"

    //echo true

    input:
        tuple val(sample_id), file(seqz_gz) from sequenza_seqz.filter{ it[1].countLines() > 1 }
        // .filter{} -- seqz files that contain more than just a header line

    output:
        tuple val(sample_id), file("${sample_id}.binned.seqz.gz") into binned_seqz

    //publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    sequenza-utils seqz_binning -w 50 -s ${seqz_gz} -o ${sample_id}.binned.seqz.gz
    """
}

process sequenza_R {
    label 'sequenza'

    tag "${sample_id}"

    //echo true

    input:
        tuple val(sample_id), file(binned_seqz_gz) from binned_seqz

    output:
        tuple val(sample_id), file("${sample_id}.nitz.cellularity.txt"), file("${sample_id}.nitz.ploidy.txt"), file("${sample_id}.nitz.copynumber_calls.txt") into sequenza_R_files
        file("${sample_id}_alternative_fit.pdf")
        file("${sample_id}_alternative_solutions.txt")
        file("${sample_id}_chromosome_depths.pdf")
        file("${sample_id}_chromosome_view.pdf")
        file("${sample_id}_CN_bars.pdf")
        file("${sample_id}.nitz.cellularity.jpg")
        file("${sample_id}_model_fit.pdf")
        file("${sample_id}_confints_CP.txt")
        file("${sample_id}_CP_contours.pdf")
        file("${sample_id}_gc_plots.pdf")
        file("${sample_id}.nitz.ave_depth.txt")
        file("${sample_id}_genome_view.pdf")

    publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
	#execute the R script
	sample_sequenza.r ${sample_id} ${binned_seqz_gz}
    """
}

process loh_score {
    label 'sequenza'

    tag "${sample_id}"

    //echo true

    input:
        path centromere_file
        tuple val(sample_id), file(cellularity), file(ploidy), file(copynumber_calls) from sequenza_R_files

    output:
        tuple val(sample_id), file("${sample_id}.nitz.score.txt") into scoring_output

    publishDir params.output, mode: 'copy', overwrite: true

    script:
    """
    LOH_score_chr_arms_V4.py ${centromere_file} ${copynumber_calls} ${sample_id}.nitz.score.txt 0.75
    echo "" >> ${sample_id}.nitz.score.txt
    echo -n "Estimated tumor cellularity: " >> ${sample_id}.nitz.score.txt
    cat ${cellularity} >> ${sample_id}.nitz.score.txt
    echo -n "Estimated ploidy: " >> ${sample_id}.nitz.score.txt
    cat ${ploidy} >> ${sample_id}.nitz.score.txt
    """
}
