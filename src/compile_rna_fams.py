#! /usr/bin/env python

import utils
from Bio import SeqIO
from Bio import AlignIO
import numpy
import seaborn
from matplotlib import pyplot
import os
import multiprocessing
import time
import datetime
import shutil

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  rfam_seed_sta_file_path = asset_dir_path + "/rfam_seed_stas_v14.3.sth"
  rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams"
  rna_seq_dir_path_4_micro_bench = asset_dir_path + "/compiled_rna_fams_4_micro_bench"
  ref_ss_dir_path = asset_dir_path + "/ref_sss"
  ref_ss_dir_path_4_micro_bench = asset_dir_path + "/ref_sss_4_micro_bench"
  if not os.path.isdir(rna_seq_dir_path):
    os.mkdir(rna_seq_dir_path)
  if not os.path.isdir(rna_seq_dir_path_4_micro_bench):
    os.mkdir(rna_seq_dir_path_4_micro_bench)
  if not os.path.isdir(ref_ss_dir_path):
    os.mkdir(ref_ss_dir_path)
  if not os.path.isdir(ref_ss_dir_path_4_micro_bench):
    os.mkdir(ref_ss_dir_path_4_micro_bench)
  max_sa_len = 200
  max_seq_num = 10
  stas = [sta for sta in AlignIO.parse(rfam_seed_sta_file_path, "stockholm") if len(sta) <= max_seq_num and len(sta[0]) <= max_sa_len and is_valid(sta)]
  num_of_stas = len(stas)
  print("# RNA families: %d" % num_of_stas)
  sample_rate = 0.02
  num_of_samples = int(sample_rate * num_of_stas)
  print("# RNA families for micro benchmark: %d" % num_of_samples)
  sampled_stas = numpy.random.choice(stas, num_of_samples, replace = False)
  for i, sta in enumerate(stas):
    css = convert_css(sta.column_annotations["secondary_structure"])
    rna_seq_file_path = os.path.join(rna_seq_dir_path, "rna_fam_%d.fa" % i)
    ref_ss_file_path = os.path.join(ref_ss_dir_path, "rna_fam_%d.fa" % i)
    rna_seq_file = open(rna_seq_file_path, "w")
    ref_ss_file = open(ref_ss_file_path, "w")
    for j, rec in enumerate(sta):
      seq_with_gaps = str(rec.seq)
      recovered_ss = recover_ss(css, seq_with_gaps)
      seq = seq_with_gaps.replace("-", "")
      rna_seq_file.write(">%d(%s)\n%s\n" % (j, rec.id, seq))
      ref_ss_file.write(">%d(%s)\n%s\n" % (j, rec.id, recovered_ss))
  for i, sta in enumerate(sampled_stas):
    css = convert_css(sta.column_annotations["secondary_structure"])
    rna_seq_file_path = os.path.join(rna_seq_dir_path_4_micro_bench, "rna_fam_%d.fa" % i)
    ref_ss_file_path = os.path.join(ref_ss_dir_path_4_micro_bench, "rna_fam_%d.fa" % i)
    rna_seq_file = open(rna_seq_file_path, "w")
    ref_ss_file = open(ref_ss_file_path, "w")
    for j, rec in enumerate(sta):
      seq_with_gaps = str(rec.seq)
      recovered_ss = recover_ss(css, seq_with_gaps)
      seq = seq_with_gaps.replace("-", "")
      rna_seq_file.write(">%d(%s)\n%s\n" % (j, rec.id, seq))
      ref_ss_file.write(">%d(%s)\n%s\n" % (j, rec.id, recovered_ss))

def is_valid(sta):
  for row in sta:
    if any(char in str(row.seq) for char in "RYWSMKHBVDN"):
      return False
  return True

def convert_css(css):
  converted_css = ""
  for char in css:
    if char == "(" or char == "<" or char == "[" or char == "{":
      converted_css += "("
    elif char == ")" or char == ">" or char == "]" or char == "}":
      converted_css += ")"
    else:
      converted_css += "."
  return converted_css

def recover_ss(css, seq_with_gaps):
  pos_map = {}
  pos = 0
  for (i, char) in enumerate(seq_with_gaps):
    if char != "-":
      pos_map[i] = pos
      pos += 1
  recovered_ss = "." * pos
  stack = []
  for (i, char) in enumerate(css):
    if char == "(":
      stack.append(i)
    elif char == ")":
      j = stack.pop()
      if seq_with_gaps[j] == "-" or seq_with_gaps[i] == "-":
        continue
      mapped_j = pos_map[j]
      mapped_i = pos_map[i]
      recovered_ss = recovered_ss[: mapped_j] + "(" + recovered_ss[mapped_j + 1 :]
      recovered_ss = recovered_ss[: mapped_i] + ")" + recovered_ss[mapped_i + 1 :]
  return recovered_ss

if __name__ == "__main__":
  main()
