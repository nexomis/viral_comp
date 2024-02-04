#!/bin/env python3

"""This is a Nextflow-compatible Python script to process FASTA files. Given a 
specific sample name, it creates separate output FASTA files for every unique 
sequence identifier while prepending a custom prefix based on the specified 
sample name. It also validates the naming convention for each sequence 
identifier ensuring only alpha-numeric characters, underscores (_), and hyphens 
(-) are used."""

import re

# Define global constants
SAMPLE_NAME = "${sample_name}"
NEW_SEQ_ID = ">${sample_name}"
FASTA_FILE = "${fasta_file}"

def validate_sequence_id(seq_id):
  pattern = r"^[A-Za-z0-9_\\-]+\$"
  if not re.match(pattern, seq_id):
    raise ValueError(f"Invalid character(s) found in sequence ID '{seq_id}'.")

def write_new_fasta(seq_name, sequences):
  if (len(sequences) > 1):
    with open(f"{seq_name}.{SAMPLE_NAME}.fa", "w") as out_file:
      out_file.write("\\n".join(sequences) + "\\n")
  return [NEW_SEQ_ID]

def main():
  sequences = [NEW_SEQ_ID]
  seq_name = ""
  validate_sequence_id(SAMPLE_NAME)
  with open(FASTA_FILE) as in_file:
    for line in in_file.readlines():
      if (line.startswith(">")):
        sequences = write_new_fasta(seq_name, sequences)
        seq_name = line.rstrip().lstrip(">")
        validate_sequence_id(seq_name)
      else:
        sequences.append(line.rstrip())
  sequences = write_new_fasta(seq_name, sequences)

main()
