#!/usr/bin/env nextflow

// Assay specific files
picard_bed_file = Channel.fromPath('/mnt/disk2/com/Genomes/Picard/MONC_ctDNA1.3_designed_probe_coords_180314_no_chr.Picard.bed')
bed_file = Channel.fromPath('/mnt/disk2/com/Genomes/BED_Files/MONC_ctDNA1.3_designed_probe_coords_180314_no_chr.bed')

// Setup the various inputs, defined in nexflow.config
fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true) //.println { it }

reference_fasta = Channel.fromPath("/mnt/disk2/com/Genomes/gatk-bundle/human_g1k_v37.fasta")
reference_index = Channel.fromPath("/mnt/disk2/com/Genomes/gatk-bundle/human_g1k_v37.fasta.{amb,ann,bwt,pac,sa,dict,fai}")

// Reference genome is used multiple times
reference_fasta.into { bwa_ref; bwa_realign_ref; picard_ref; qc_ref; filter_con_ref ; vardict_ref}
reference_index.into { bwa_ref_index; bwa_realign_ref_index; picard_ref_index; qc_ref_index; filter_con_ref_index; vardict_ref_index }

//  memory "32GB"

process bwa {
  // Align fastqs
  label 'bwa'
  tag "${sample_id}"
  input:
    file(reference_fasta) from bwa_ref
    file("*") from bwa_ref_index.collect()
    set sample_id, file(fastq1), file(fastq2) from fastq_pair_ch

  output:
    set val(sample_id), file('*.sam') into align_ch

  publishDir params.output, overwrite: true

  cpus 8
  
  script:
  """ 
  bwa mem \
  -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
  -K 10000000 \
  -C \
  -Y \
  -t ${task.cpus} \
  ${reference_fasta} \
  ${fastq1} ${fastq2} > ${sample_id}.sam
  """
}
