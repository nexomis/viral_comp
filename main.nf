#!/usr/bin/env nextflow

if ( params.help ) {
  help = """
    |viral_comp: A workflow for the comparison of viral sequences.
    |
    |Required arguments:
    |  --sample_sheet   Sample sheet file (CSV format ,) with header as below:
    |                   `sample_name,fasta_file`
    |  --ref_sample     Reference sample for the annotation
    |  --out_dir        Output directory
    |
  """.stripMargin()
  // Print the help with the stripped margin and exit
  println(help)
  exit(0)
}

// Validate input parameters
if (params.sample_sheet == null) {
  error """No sample sheet provided. Use --sample_sheet to specify the path to the 
          | sample sheet file.""".stripMargin()
}
if (params.ref_sample == null) {
  error """No ref_sample provided. Use --ref_sample to specify the reference
          | sample for the annotation.""".stripMargin()
}

if (params.out_dir == null) {
  error "No out_dir provided. Use --out_dir to specify the output directory"
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
  tuple val(seq_name), path("${seq_name}.fa")

  script:
  """
  #!/bin/bash
  cat *.fa > ${seq_name}.fa
  """

}

process ANNOT_REF {

  container "staphb/prokka:latest"

  label 'multi_cpu'
  label 'low_dyn_mem'

  input:
  tuple val(seq_name), path(fasta_files)

  output:
  tuple val(seq_name), path("reference_annotation/*")

  publishDir "$params.out_dir"

  script:
  """
  #!/bin/bash
  prokka --kingdom Viruses --cpu $task.cpus --outdir reference_annotation \\
  --locustag L --cdsrnaolap --quiet $fasta_files
  """

}

process ALIGN_SEGMENT {
  container "cbcrg/tcoffee@sha256:9b42d8ca64720e44a12ca2d56cad7a897b3dfa9ebd1cefae2b651d343d896461"

  label 'single_cpu'
  label 'low_mem'

  input:
  tuple val(seq_name), path(fasta_file)

  output:
  tuple val(seq_name), path("${seq_name}.aln.fa")

  script:
  """
  #!/bin/bash
  t_coffee $fasta_file -type=dna -mode mcoffee -outfile ${seq_name}.aln.fa -thread $task.cpus -debug -output fasta_aln
  """
}

process DISTRIBUTE_SEQS {
  container 'python:3.10'

  label 'single_cpu'
  label 'low_mem'

  input:
  tuple val(ref_name), path(annot_file), val(seq_name), path(aln_file)

  output:
  path "{cds,nc}/*.fa", optional: true

  publishDir "$params.out_dir", pattern: "nc/*"

  script:
  template "parse_first_aln_with_annot.py"

}

process ALIGN_CDS {
  container 'ghcr.io/nexomis/macse:v2.07.0'

  label 'single_cpu'
  label 'low_mem'

  input:
  tuple val(seq_name), path(fasta_file)

  output:
  tuple val(seq_name), path("aligned_cds/${seq_name}/NT_aln.fa"), path("aligned_cds/${seq_name}/AA_aln.fa")

  publishDir "$params.out_dir"

  script:
  """
  #!/bin/bash
  mkdir -p aligned_cds
  mkdir -p aligned_cds/${seq_name}
  macse -prog alignSequences -seq $fasta_file -out_NT aligned_cds/${seq_name}/NT_aln.raw.fa -out_AA aligned_cds/${seq_name}/AA_aln.raw.fa
  fold -w 80 aligned_cds/${seq_name}/NT_aln.raw.fa > aligned_cds/${seq_name}/NT_aln.fa
  fold -w 80 aligned_cds/${seq_name}/AA_aln.raw.fa > aligned_cds/${seq_name}/AA_aln.fa
  rm aligned_cds/${seq_name}/NT_aln.raw.fa aligned_cds/${seq_name}/AA_aln.raw.fa
  """
}

process FORMAT_ALN_CDS {
  container 'python:3.10'

  label 'single_cpu'
  label 'low_mem'

  input:
  tuple val(ref_sample), val(seq_name), path(nt_aln_file), path(aa_aln_file)

  output:
  path("cds_mutable/*/*"), optional: true

  publishDir "$params.out_dir"

  script:
  template "format_aligned_cds.py"
}

process FORMAT_ALN_NC {
  container 'python:3.10'

  label 'single_cpu'
  label 'low_mem'

  input:
  tuple val(ref_sample), val(seq_name), path(nt_aln_file)

  output:
  path("nc_mutable/*"), optional: true

  publishDir "$params.out_dir"

  script:
  template "format_aligned_nc.py"
}

workflow{

  Channel.fromPath(params.sample_sheet)
  | splitCsv(header: true, sep: ',', strip: true)
  | map { row -> tuple(row.sample_name, row.fasta_file) }
  | filter { it[0] == params.ref_sample }
  | set {refInput}

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
  | ALIGN_SEGMENT
  | set {concatFasta}

  ANNOT_REF(refInput)
  | flatMap { key, items ->
        items.collect { item -> [key, item] }
    }
  | filter { it[1].getExtension() == "gff" }
  | combine(concatFasta)
  | set {seqsToDistribute}

  DISTRIBUTE_SEQS(seqsToDistribute)
  | flatten()
  | map {tuple( it.parent.name, it.name.split('\\.')[0], it )}
  | branch {
      cds: it[0] == "cds"
      nc: it[0] == "nc"
    }
  | set { ditributedSeqs }

  ditributedSeqs.cds 
  | map {tuple(it[1], it[2])}
  | ALIGN_CDS
  | map {tuple(params.ref_sample, it[0], it[1], it[2])}
  | FORMAT_ALN_CDS

  ditributedSeqs.nc
  | map {tuple(params.ref_sample, it[1], it[2])}
  | FORMAT_ALN_NC

}



