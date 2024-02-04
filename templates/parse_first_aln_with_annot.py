#!/bin/env python3 

import re
import os
import sys

SEQ_NAME = "${seq_name}"
ALN_FILE = "${aln_file}"
ANNOT_FILE = "${annot_file}"
REF_NAME = "${ref_name}"
COMP = { "A":"T", "C":"G", "G":"C", "T":"A", "U":"A", "M":"K", "R":"Y", "W":"W", 
"S":"S", "Y":"R", "K":"M", "V":"B", "H":"D", "D":"H", "B":"V", "N":"N", 
"-":"-", ".":".", "!":"!", "?":"?" }

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

def parse_annot(seq_name, annot_file):
  i = 1
  annot = []
  existing_names = []
  with open(annot_file, "r") as file:
    for line in file.readlines():
      if line.startswith(">"): break
      if line.startswith("#"): continue
      sline = line.split("\\t")
      line_seq_name = sline[0]
      type_seq = sline[2]
      if type_seq != "CDS": continue

      # Extract name if CDS
      match = re.search(r".*;Name=([^;]+);.*", sline[8])
      if match:
        cds_name = match.groups()[0]
      else:
        cds_name = line_seq_name + "_" + str(i)
        i+=1

      # Make sure name is unique (works if no splicing)
      j = 1
      while cds_name in existing_names:
        new_name = cds_name + "_" + str(j)
        if new_name in existing_names:
          j+=1
        else: 
          cds_name = new_name
      existing_names.append(cds_name)

      # Focus only in current seq
      if line_seq_name != seq_name: continue
      annot.append((int(sline[3]), int(sline[4]), sline[6], cds_name))
  return sorted(annot, key=lambda x: x[0])

def make_ref2aln_dict(ref_name, seq_dict):
  ref2aln=[0]
  i=0
  for char in seq_dict[ref_name]:
    i+=1
    if char != "-":
      ref2aln.append(i)
  return ref2aln

def get_non_coding(annot, max_len):
  starts = [row[0] for row in annot]
  ends = [row[1] for row in annot]
  names = [row[3] for row in annot]
  starts.append(max_len + 1)
  ends.append(max_len + 1)
  last_end = 0
  last_name = SEQ_NAME + "_5p"
  names.append(SEQ_NAME + "_3p")
  nc_annot = []
  for i in range(len(starts)):
    nc_annot.append([
      last_end + 1, starts[i]-1, "+", last_name + "_to_" + names[i]
    ])
    last_end = ends[i]
    last_name = names[i]
  return nc_annot

def tr_annot(annot, ref2aln):
  tr_annot = []
  for cds in annot:
    tr_annot.append((ref2aln[cds[0]], ref2aln[cds[1]] ,cds[2] ,cds[3]))
  return tr_annot

def get_seqs_from_annot(annot, seq_dict, remove_gap):
  annot_seq_dict = dict()
  for cds in annot:
    sub_seq_dict = dict()
    for seq_name, seq in seq_dict.items():
      sub_seq = seq[cds[0]-1:cds[1]].upper()
      if remove_gap: sub_seq = re.sub("[!?.-]", "", sub_seq).strip()
      if (cds[2] == "-"):
        sub_seq = "".join([COMP[base] for base in reversed(sub_seq)])
      sub_seq_dict[seq_name] = sub_seq
    annot_seq_dict[cds[3]] = sub_seq_dict
  return annot_seq_dict

def write_seqs(out_dir, annot_seq_dict):
  os.makedirs(out_dir, exist_ok = True)
  for nc, seq_dict in annot_seq_dict.items():
    filename = f"{out_dir}/{nc}.fa"
    with open(filename, 'w') as fh:
      for seq_name, seq in seq_dict.items():
        fh.write(f">{seq_name}\\n")
        for i in range(0, len(seq), 80):
          fh.write(seq[i:i+80].upper()+"\\n")

def main():
  seq_dict = parse_aln(ALN_FILE)
  annot = parse_annot(SEQ_NAME, ANNOT_FILE)
  ref2aln = make_ref2aln_dict(REF_NAME, seq_dict)
  annot = tr_annot(annot, ref2aln)
  nc_annot = get_non_coding(annot, len(seq_dict[REF_NAME]))
  annot_seq_dict = get_seqs_from_annot(annot, seq_dict, True)
  nc_annot_seq_dict = get_seqs_from_annot(nc_annot, seq_dict, False)
  write_seqs("cds", annot_seq_dict)
  write_seqs("nc", nc_annot_seq_dict)

main()
