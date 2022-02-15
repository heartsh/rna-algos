#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os
import multiprocessing
import time
import datetime
import shutil
from os import path

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  temp_dir_path = "/tmp/run_ss_estimation_programs_%s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
  if not os.path.isdir(temp_dir_path):
    os.mkdir(temp_dir_path)
  gammas = [2. ** i for i in range(-7, 11)]
  centroid_estimator_turner_params = []
  centroid_estimator_contra_params = []
  centroidfold_params_turner = []
  centroidfold_params_contra = []
  centroid_estimator_turner_dir_path = asset_dir_path + "/centroid_estimator_turner"
  centroid_estimator_contra_dir_path = asset_dir_path + "/centroid_estimator_contra"
  centroidfold_dir_path_turner = asset_dir_path + "/centroidfold_turner"
  centroidfold_dir_path_contra = asset_dir_path + "/centroidfold_contra"
  if not os.path.isdir(centroid_estimator_turner_dir_path):
    os.mkdir(centroid_estimator_turner_dir_path)
  if not os.path.isdir(centroid_estimator_contra_dir_path):
    os.mkdir(centroid_estimator_contra_dir_path)
  if not os.path.isdir(centroidfold_dir_path_turner):
    os.mkdir(centroidfold_dir_path_turner)
  if not os.path.isdir(centroidfold_dir_path_contra):
    os.mkdir(centroidfold_dir_path_contra)
  rna_dir_path = asset_dir_path + "/compiled_rna_fams"
  # rna_dir_path = asset_dir_path + "/compiled_rna_fams_4_micro_bench"
  sub_thread_num = 4
  for rna_file in os.listdir(rna_dir_path):
    if not rna_file.endswith(".fa"):
      continue
    rna_file_path = os.path.join(rna_dir_path, rna_file)
    (rna_family_name, extension) = os.path.splitext(rna_file)
    centroid_estimator_turner_output_dir_path = os.path.join(centroid_estimator_turner_dir_path, rna_family_name)
    centroid_estimator_contra_output_dir_path = os.path.join(centroid_estimator_contra_dir_path, rna_family_name)
    centroidfold_output_dir_path_turner = os.path.join(centroidfold_dir_path_turner, rna_family_name)
    centroidfold_output_dir_path_contra = os.path.join(centroidfold_dir_path_contra, rna_family_name)
    if not os.path.isdir(centroid_estimator_turner_output_dir_path):
      os.mkdir(centroid_estimator_turner_output_dir_path)
    if not os.path.isdir(centroid_estimator_contra_output_dir_path):
      os.mkdir(centroid_estimator_contra_output_dir_path)
    if not os.path.isdir(centroidfold_output_dir_path_turner):
      os.mkdir(centroidfold_output_dir_path_turner)
    if not os.path.isdir(centroidfold_output_dir_path_contra):
      os.mkdir(centroidfold_output_dir_path_contra)
    centroid_estimator_turner_command = "centroid_estimator -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + centroid_estimator_turner_output_dir_path
    centroid_estimator_turner_params.insert(0, centroid_estimator_turner_command)
    centroid_estimator_contra_command = "centroid_estimator -c -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + centroid_estimator_contra_output_dir_path
    centroid_estimator_contra_params.insert(0, centroid_estimator_contra_command)
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      output_file = "gamma=" + gamma_str + ".fa"
      centroidfold_output_file_path_turner = os.path.join(centroidfold_output_dir_path_turner, output_file)
      centroidfold_params_turner.insert(0, (rna_file_path, centroidfold_output_file_path_turner, gamma_str, False))
      centroidfold_output_file_path_contra = os.path.join(centroidfold_output_dir_path_contra, output_file)
      centroidfold_params_contra.insert(0, (rna_file_path, centroidfold_output_file_path_contra, gamma_str, True))
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  begin = time.time()
  pool.map(utils.run_command, centroid_estimator_turner_params)
  centroid_estimator_turner_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(utils.run_command, centroid_estimator_contra_params)
  centroid_estimator_contra_elapsed_time = time.time() - begin
  pool = multiprocessing.Pool(num_of_threads)
  begin = time.time()
  pool.map(run_centroidfold, centroidfold_params_turner)
  centroidfold_elapsed_time_turner = time.time() - begin
  begin = time.time()
  pool.map(run_centroidfold, centroidfold_params_contra)
  centroidfold_elapsed_time_contra = time.time() - begin
  print("The elapsed time of centroid estimator (Turner) = %f [s]." % centroid_estimator_turner_elapsed_time)
  print("The elapsed time of centroid estimator (CONTRAfold) = %f [s]." % centroid_estimator_contra_elapsed_time)
  print("The elapsed time of CentroidFold (Turner) = %f [s]." % centroidfold_elapsed_time_turner)
  print("The elapsed time of CentroidFold (CONTRAfold) = %f [s]." % centroidfold_elapsed_time_contra)
  shutil.rmtree(temp_dir_path)

def run_centroidfold(centroidfold_params):
  (rna_file_path, centroidfold_output_file_path, gamma_str, is_contra) = centroidfold_params
  centroidfold_command = "centroid_fold --engine %s " % ("CONTRAfold" if is_contra else "McCaskill") + rna_file_path + " -g " + gamma_str
  (output, _, _) = utils.run_command(centroidfold_command)
  lines = [line.split()[0] for (i, line) in enumerate(str(output).split("\\n")) if i % 3 == 2]
  centroidfold_output_file = open(centroidfold_output_file_path, "w+")
  centroidfold_output_buf = ""
  for (i, line) in enumerate(lines):
    centroidfold_output_buf += ">%d\n%s\n\n" % (i, line)
  centroidfold_output_file.write(centroidfold_output_buf)
  centroidfold_output_file.close()

if __name__ == "__main__":
  main()
