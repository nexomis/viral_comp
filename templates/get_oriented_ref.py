#!/bin/env python3 

import re
import os
import sys

SEQ_NAME = "${seq_name}"
ALN_FILE = "${aln_file}"
REF_SAMPLE = "${ref_sample}"

def parse_aln(aln_file):
  seq_dict = dict()
  seq = ""
  sample_name = ""
  with open(aln_file, "r") as file:
    for line in file.readlines():
      if (line.startswith(">")):
        if sample_name != "":
          seq_dict[sample_name] = seq
        seq = ""
        sample_name = line.rstrip().lstrip(">")
      else:
        seq += line.rstrip() 
    if sample_name != "":
      seq_dict[sample_name] = seq
  return seq_dict

def main():
  seq_dict = parse_aln(ALN_FILE)
  ref_seq = seq_dict[REF_SAMPLE].replace("-", "")
  with open(f"oriented_{REF_SAMPLE}_{SEQ_NAME}.fa", "w") as fh:
    fh.write(f">{SEQ_NAME}\\n")
    for i in range(0, len(ref_seq), 80):
      fh.write(ref_seq[i:i+80].upper()+"\\n")

main()
