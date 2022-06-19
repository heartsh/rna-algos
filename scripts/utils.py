import os
import matplotlib
from matplotlib import pylab
import numpy
import subprocess
from Bio import SeqIO

bracket_pairs = [("(", ")"), ("<", ">"), ("{", "}"), ("[", "]"), ("A", "a"), ("B", "b"), ("C", "c"), ("D", "d"), ("E", "e"), ]

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)

def get_sss(ss_file_path):
  sss = []
  ss_strings = []
  ss_strings = [rec.seq for rec in SeqIO.parse(ss_file_path, "fasta")]
  sss = []
  for (i, ss_string) in enumerate(ss_strings):
    sss.append({})
    for (left, right) in bracket_pairs:
      stack = []
      for (j, char) in enumerate(ss_string):
        if char == left:
          stack.append(j)
        elif char == right:
          pos = stack.pop()
          sss[i][(pos, j)] = True
  return sss

def run_command(command):
  subproc = subprocess.Popen(
    command,
    stdout = subprocess.PIPE,
    shell = True
  )
  (output, error) = subproc.communicate()
  returned_code = subproc.wait()
  return (output, error, returned_code)
