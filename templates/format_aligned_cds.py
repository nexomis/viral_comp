#!/bin/env python3

import os
import sys

SEQ_NAME = "${seq_name}" 
# Alignement with amino acid
AA_FASTA_FILE = "${aa_aln_file}"
# Same alignement with nucleotid
NT_FASTA_FILE = "${nt_aln_file}"
# reference sample
REF_SAMPLE = "${ref_sample}"

STOP_CODONS = ["TAG", "TAA", "TGA"]

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

def compare_sequences(aa_seqs, nt_seqs, ref_sample):
  ref_aa_seq = aa_seqs[ref_sample]
  ref_nt_seq = nt_seqs[ref_sample]
  aa_pos_ref = -1
  data_pos_all_aa = [["ref_aa_pos", "aln_aa_pos"]]
  data_pos_all_nt = [["ref_nt_pos", "aln_nt_pos"]]
  data_pos_ns_aa = [["ref_aa_pos", "aln_aa_pos"]]
  data_pos_ns_nt = [["ref_nt_pos", "aln_nt_pos"]]
  all_samples = [ref_sample]
  for sample_name in sorted(aa_seqs.keys()):
    if sample_name != ref_sample:
      all_samples.append(sample_name)
  data_spl_all_nt = [[ spl + "_nt" for spl in all_samples]]
  data_spl_ns_nt = [[ spl + "_nt" for spl in all_samples]]
  data_spl_all_aa = [[ spl + "_aa" for spl in all_samples]]
  data_spl_ns_aa = [[ spl + "_aa" for spl in all_samples]]
  data_spl_all_mut = [[ spl + "_mut" for spl in all_samples]]
  data_spl_ns_mut = [[ spl + "_mut" for spl in all_samples]]
  for aa_pos in range(len(ref_aa_seq)):
    codon_pos = aa_pos * 3
    ref_aa = ref_aa_seq[aa_pos]
    if not ref_aa in ["-", "?", "!", "*"]:
      aa_pos_ref += 1
    ref_codon = ref_nt_seq[codon_pos:codon_pos + 3]
    has_nt_diff = False
    has_aa_diff = False
    codons = [ref_codon]
    muts = ["."]
    aas = [ref_aa]
    for sample_name in sorted(aa_seqs.keys()):
      if sample_name != ref_sample:
        codon = nt_seqs[sample_name][codon_pos:codon_pos + 3]
        aa = aa_seqs[sample_name][aa_pos]
        if codon == ref_codon: 
          codons.append(".")
          muts.append(".")
          aas.append(".")
        else :
          has_nt_diff = True
          codons.append(codon)
          if aa == ref_aa:
            aas.append(".")
            muts.append("SYN")
          else:
            has_aa_diff = True
            aas.append(aa)
            if ref_aa == "-":
              muts.append("INS")
            elif aa in ["!", "?"] or ref_aa in ["!", "?"]:
              muts.append("FS")
            elif aa == "-":
              muts.append("DEL")
            elif aa in STOP_CODONS and (not ref_aa in STOP_CODONS):
              muts.append("STOP")
            elif ref_aa in STOP_CODONS and (not aa in STOP_CODONS):
              muts.append("STOP_LOSS")
            else:
              muts.append("MIS")
    if has_nt_diff:
      data_pos_all_aa.append([str(aa_pos_ref+1), str(aa_pos+1)])
      data_pos_all_nt.append([str(aa_pos_ref*3+3), str(aa_pos*3+3)])
      data_spl_all_nt.append(codons)
      data_spl_all_aa.append(aas)
      data_spl_all_mut.append(muts)
    if has_aa_diff:
      data_pos_ns_aa.append([str(aa_pos_ref+1), str(aa_pos_ref+1)])
      data_pos_ns_nt.append([str(aa_pos_ref*3+3), str(aa_pos_ref*3+3)])
      data_spl_ns_nt.append(codons)
      data_spl_ns_aa.append(aas)
      data_spl_ns_mut.append(muts)
  os.makedirs("cds_mutable", exist_ok = True)
  if (len(data_pos_all_aa)) > 1:
    os.makedirs("cds_mutable/all_muts", exist_ok = True)
    with open(f"cds_mutable/all_muts/{SEQ_NAME}.tsv", "w") as fh:
      for i in range(len(data_pos_all_aa)):
        fh.write("\\t".join(data_pos_all_aa[i]) + "\\t")
        fh.write("\\t".join(data_pos_all_nt[i]) + "\\t")
        fh.write("\\t".join(data_spl_all_mut[i]) + "\\t")
        fh.write("\\t".join(data_spl_all_nt[i]) + "\\t")
        fh.write("\\t".join(data_spl_all_aa[i]) + "\\n")
  if (len(data_pos_ns_aa)) > 1:
    os.makedirs("cds_mutable/ns_muts", exist_ok = True)
    os.makedirs("cds_mutable/ns_muts_short", exist_ok = True)
    with open(f"cds_mutable/ns_muts/{SEQ_NAME}.tsv", "w") as fh:
      for i in range(len(data_pos_ns_aa)):
        fh.write("\\t".join(data_pos_ns_aa[i]) + "\\t")
        fh.write("\\t".join(data_pos_ns_nt[i]) + "\\t")
        fh.write("\\t".join(data_spl_ns_mut[i]) + "\\t")
        fh.write("\\t".join(data_spl_ns_nt[i]) + "\\t")
        fh.write("\\t".join(data_spl_ns_aa[i]) + "\\n")
    with open(f"cds_mutable/ns_muts_short/{SEQ_NAME}.tsv", "w") as fh:
      data_spl_ns_aa[0] = [x.rstrip("_aa") for x in data_spl_ns_aa[0]]
      for i in range(len(data_pos_ns_aa)):
        fh.write("\\t".join(data_pos_ns_aa[i]) + "\\t")
        fh.write("\\t".join(data_spl_ns_aa[i]) + "\\n")

def main():
  aa_seqs = parse_aln(AA_FASTA_FILE)
  nt_seqs = parse_aln(NT_FASTA_FILE)
  compare_sequences(aa_seqs, nt_seqs, REF_SAMPLE)

main()
