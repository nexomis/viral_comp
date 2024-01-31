#!/bin/env python3

SAMPLE_NAME = "${sample_name}"
NEW_SEQ_ID = ">${sample_name}"
FASTA_FILE = ${fasta_file}

def write_new_fasta(seq_name, sequences):
  if (len(sequences) > 1):
    with open(seq_name + "." + SAMPLE_NAME + ".fa", "w") as out_file:
      out_file.writelines(sequences)
  return([NEW_SEQ_ID])

def main():
  sequences = [NEW_SEQ_ID]
  with open(FASTA_FILE) as in_file:
    for line in in_file.readlines():
      if (line.startswith(">")):
        sequences = write_new_fasta(seq_name, sequences)
        seq_name = line.rstrip().lstrip(">")
      else:
        sequences.append(line.rstrip())
  sequences= write_new_fasta(seq_name, sequences)

