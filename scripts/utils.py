import os
import matplotlib
from matplotlib import pylab
import numpy
import subprocess
from Bio import SeqIO

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)

def get_ss_and_flat_ss(ss_string):
  ss = {}
  flat_ss = {}
  stack = []
  for (i, char) in enumerate(ss_string):
    if char == "(":
      stack.append(i)
    elif char == ")":
      pos = stack.pop()
      ss[(pos, i)] = True
      flat_ss[pos] = True
      flat_ss[i] = True
  return ss, flat_ss

def get_ss_strings(ss_file_path):
  ss_strings = [filter(rec.seq) for rec in SeqIO.parse(ss_file_path, "fasta")]
  return ss_strings

def filter(seq):
  new_seq = str(seq).replace("A", "").replace("C", "").replace("G", "").replace("U", "")
  return new_seq

def get_sss_and_flat_sss(ss_strings):
  return list(map(get_ss_and_flat_ss, ss_strings))

def run_command(command):
  subproc = subprocess.Popen(
    command,
    stdout = subprocess.PIPE,
    shell = True
  )
  (output, error) = subproc.communicate()
  returned_code = subproc.wait()
  return (output, error, returned_code)

def get_bpp_mats(bpp_mat_file_path, seq_lens):
  bpp_mats = {}
  bpp_mat_file = open(bpp_mat_file_path)
  lines = bpp_mat_file.readlines()
  lines = [line for line in lines if line[0].isdigit() or line[0].startswith(">")]
  num_of_lines = len(lines)
  for i in range(0, num_of_lines - 1, 2):
    rna_id = int(lines[i][1 :])
    seq_len = seq_lens[rna_id]
    bpp_mat = numpy.zeros((seq_len, seq_len))
    for string in lines[i + 1].strip().split(" "):
      substrings = string.split(",")
      (j, k, bpp) = (int(substrings[0]), int(substrings[1]), min(1, float(substrings[2])))
      bpp_mat[j, k] = bpp
    bpp_mats[rna_id] = bpp_mat
  return bpp_mats

def get_upp_mats(upp_mat_file_path, seq_lens):
  upp_mats = {}
  upp_mat_file = open(upp_mat_file_path)
  lines = upp_mat_file.readlines()
  lines = [line for line in lines if line[0].isdigit() or line[0].startswith(">")]
  num_of_lines = len(lines)
  for i in range(0, num_of_lines - 1, 2):
    rna_id = int(lines[i][1 :])
    seq_len = seq_lens[rna_id]
    upp_mat = numpy.zeros(seq_len)
    for string in lines[i + 1].strip().split(" "):
      substrings = string.split(",")
      (j, upp) = (int(substrings[0]), min(1, float(substrings[1])))
      upp_mat[j] = upp
    upp_mats[rna_id] = upp_mat
  return upp_mats

def get_capr_prof_seqs(capr_output_file_path):
  capr_prof_seqs = {}
  capr_output_file = open(capr_output_file_path)
  lines = capr_output_file.readlines()
  lines = [line for line in lines]
  num_of_lines = len(lines)
  for line in lines[1:]:
    strings = line.split()
    prof_type = strings[0]
    profs = []
    for string in strings[1:]:
      profs.append(float(string))
    capr_prof_seqs[prof_type] = profs
  return capr_prof_seqs
