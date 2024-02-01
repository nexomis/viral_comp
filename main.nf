#!/usr/bin/env nextflow

if ( params.help ) {
  help = """
    |viral_comp: A workflow for the comparison of viral sequences.
    |
    |Required arguments:
    |  --sample_sheet   Sample sheet file (CSV format ,) with header as below:
    |                   `sample_name,fasta_file`
    |  --ref_sample     Reference sample for the annotation
    |
  """.stripMargin()
  // Print the help with the stripped margin and exit
  println(help)
  exit(0)
}

// Validate input parameters
if (params.sample_sheet == null) {
  error """No samplesheet provided. Use --samplesheet to specify the path to the 
          | samplesheet file.""".stripMargin()
}
if (params.ref_sample == null) {
  error """No ref_sample provided. Use --ref_sample reference sample for the
          | annotation.""".stripMargin()
}

// Process 1: Split each fasta file into separate files per sequence
process SPLIT_FASTA {
  container 'python:3.10'

  label 'single_cpu'
  label 'low_mem'

  input:
  tuple val(sample_name), path(fasta_file)

  output:
  path("*.fa")

  script:
  template "split_fasta_per_sequence.py"

}

process CONCAT_FASTA {

  container "debian:stable-slim"

  label 'single_cpu'
  label 'low_mem'

  input:
  tuple val(seq_name), path(fasta_files)

  output:
  path "${seq_name}.fa"

  script:
  """
  #!/bin/bash
  cat *.fa > ${seq_name}.fa
  """

}

workflow{

  Channel.fromPath(params.sample_sheet)
  | splitCsv(header: true, sep: ',', strip: true)
  | map { row -> tuple(row.sample_name, row.fasta_file) }
  | set { samplesInput }

  SPLIT_FASTA(samplesInput)
  | flatten() 
  | map {tuple( it.name.split('\\.')[0], it )}
  | groupTuple(by: 0)
  | set {splitFastas}

  CONCAT_FASTA(splitFastas)
  | view()


}

