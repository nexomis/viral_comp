#!/bin/env python3

import os
import sys

SEQ_NAME = "${seq_name}" 
# Same alignement with nucleotid
NT_FASTA_FILE = "${nt_aln_file}"
# reference sample
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
      seq_dict[sample_name] = seq.upper()
  return seq_dict

def compare_sequences(nt_seqs, ref_sample):
  ref_nt_seq = nt_seqs[ref_sample]
  nt_pos_ref = -1
  data_pos_nt = [["ref_nt_pos", "aln_nt_pos"]]
  all_samples = [ref_sample]
  for sample_name in sorted(nt_seqs.keys()):
    if sample_name != ref_sample:
      all_samples.append(sample_name)
  data_spl_nt = [all_samples]
  for nt_pos in range(len(ref_nt_seq)):
    ref_nt = ref_nt_seq[nt_pos]
    if not ref_nt in ["-", "?", "!", "*"]:
      nt_pos_ref += 1
    has_nt_diff = False
    nts = [ref_nt]
    for sample_name in sorted(nt_seqs.keys()):
      if sample_name != ref_sample:
        nt = nt_seqs[sample_name][nt_pos]
        if nt == ref_nt: 
          nts.append(".")
        else :
          has_nt_diff = True
          nts.append(nt)
    if has_nt_diff:
      data_pos_nt.append([str(nt_pos_ref+1), str(nt_pos+1)])
      data_spl_nt.append(nts)
  os.makedirs("nc_mutable", exist_ok = True)
  if (len(data_pos_nt)) > 1:
    with open(f"nc_mutable/{SEQ_NAME}.tsv", "w") as fh:
      for i in range(len(data_pos_nt)):
        fh.write("\\t".join(data_pos_nt[i]) + "\\t")
        fh.write("\\t".join(data_spl_nt[i]) + "\\n")

def main():
  nt_seqs = parse_aln(NT_FASTA_FILE)
  compare_sequences(nt_seqs, REF_SAMPLE)

main()
