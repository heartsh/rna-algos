#! /usr/bin/env python

import utils
from Bio import SeqIO
import seaborn
from matplotlib import pyplot
import os
from sklearn.metrics import roc_curve
import math
from math import sqrt
import multiprocessing

seaborn.set()

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  centroid_estimator_turner_ppvs = []
  centroid_estimator_turner_senss = []
  centroid_estimator_turner_fprs = []
  centroid_estimator_turner_f1_scores = []
  centroid_estimator_turner_mccs = []
  centroid_estimator_contra_ppvs = []
  centroid_estimator_contra_senss = []
  centroid_estimator_contra_fprs = []
  centroid_estimator_contra_f1_scores = []
  centroid_estimator_contra_mccs = []
  centroidfold_ppvs = []
  centroidfold_senss = []
  centroidfold_fprs = []
  centroidfold_f1_scores = []
  centroidfold_mccs = []
  gammas = [2. ** i for i in range(-7, 11)]
  centroid_estimator_turner_ss_dir_path = asset_dir_path + "/centroid_estimator_turner"
  centroid_estimator_contra_ss_dir_path = asset_dir_path + "/centroid_estimator_contra"
  centroidfold_ss_dir_path = asset_dir_path + "/centroidfold"
  rna_fam_dir_path = asset_dir_path + "/ref_sss"
  # rna_fam_dir_path = asset_dir_path + "/ref_sss_4_micro_bench"
  pool = multiprocessing.Pool(num_of_threads)
  for gamma in gammas:
    centroid_estimator_turner_count_params = []
    centroid_estimator_contra_count_params = []
    centroidfold_count_params = []
    gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_ss_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      ref_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(ref_ss_file_path))
      centroid_estimator_turner_estimated_ss_dir_path = os.path.join(centroid_estimator_turner_ss_dir_path, rna_fam_name)
      centroid_estimator_turner_estimated_ss_file_path = os.path.join(centroid_estimator_turner_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      centroid_estimator_turner_count_params.insert(0, (centroid_estimator_turner_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      centroid_estimator_contra_estimated_ss_dir_path = os.path.join(centroid_estimator_contra_ss_dir_path, rna_fam_name)
      centroid_estimator_contra_estimated_ss_file_path = os.path.join(centroid_estimator_contra_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      centroid_estimator_contra_count_params.insert(0, (centroid_estimator_contra_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      centroidfold_estimated_ss_dir_path = os.path.join(centroidfold_ss_dir_path, rna_fam_name)
      centroidfold_estimated_ss_file_path = os.path.join(centroidfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      centroidfold_count_params.insert(0, (centroidfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
    results = pool.map(get_pos_neg_counts, centroid_estimator_turner_count_params)
    tp, tn, fp, fn = final_sum(results)
    ppv = get_ppv(tp, fp)
    sens = get_sens(tp, fn)
    fpr = get_fpr(tn, fp)
    centroid_estimator_turner_ppvs.insert(0, ppv)
    centroid_estimator_turner_senss.insert(0, sens)
    centroid_estimator_turner_fprs.insert(0, fpr)
    centroid_estimator_turner_f1_scores.append(get_f1_score(ppv, sens))
    centroid_estimator_turner_mccs.append(get_mcc(tp, tn, fp, fn))
    results = pool.map(get_pos_neg_counts, centroid_estimator_contra_count_params)
    tp, tn, fp, fn = final_sum(results)
    ppv = get_ppv(tp, fp)
    sens = get_sens(tp, fn)
    fpr = get_fpr(tn, fp)
    centroid_estimator_contra_ppvs.insert(0, ppv)
    centroid_estimator_contra_senss.insert(0, sens)
    centroid_estimator_contra_fprs.insert(0, fpr)
    centroid_estimator_contra_f1_scores.append(get_f1_score(ppv, sens))
    centroid_estimator_contra_mccs.append(get_mcc(tp, tn, fp, fn))
    results = pool.map(get_pos_neg_counts, centroidfold_count_params)
    tp, tn, fp, fn = final_sum(results)
    ppv = get_ppv(tp, fp)
    sens = get_sens(tp, fn)
    fpr = get_fpr(tn, fp)
    centroidfold_ppvs.insert(0, ppv)
    centroidfold_senss.insert(0, sens)
    centroidfold_fprs.insert(0, fpr)
    centroidfold_f1_scores.append(get_f1_score(ppv, sens))
    centroidfold_mccs.append(get_mcc(tp, tn, fp, fn))
  line_1, = pyplot.plot(centroid_estimator_turner_ppvs, centroid_estimator_turner_senss, label = "Centroid estimator (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(centroid_estimator_contra_ppvs, centroid_estimator_contra_senss, label = "Centroid estimator (CONTRAfold)", marker = "p", linestyle = ":")
  line_3, = pyplot.plot(centroidfold_ppvs, centroidfold_senss, label = "CentroidFold (CONTRAfold)", marker = "s", linestyle = "--")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  pyplot.legend(handles = [line_1, line_2, line_3], loc = "lower left")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  line_1, = pyplot.plot(centroid_estimator_turner_fprs, centroid_estimator_turner_senss, label = "Centroid estimator (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(centroid_estimator_contra_fprs, centroid_estimator_contra_senss, label = "Centroid estimator (CONTRAfold)", marker = "p", linestyle = ":")
  line_3, = pyplot.plot(centroidfold_fprs, centroidfold_senss, label = "CentroidFold (CONTRAfold)", marker = "s", linestyle = "--")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  gammas = [i for i in range(-7, 11)]
  line_1, = pyplot.plot(gammas, centroid_estimator_turner_f1_scores, label = "Centroid estimator (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, centroid_estimator_contra_f1_scores, label = "Centroid estimator (CONTRAfold)", marker = "p", linestyle = ":")
  line_3, = pyplot.plot(gammas, centroidfold_f1_scores, label = "CentroidFold (CONTRAfold)", marker = "s", linestyle = "--")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.legend(handles = [line_1, line_2, line_3], loc = "lower right")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  line_1, = pyplot.plot(gammas, centroid_estimator_turner_mccs, label = "Centroid estimator (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, centroid_estimator_contra_mccs, label = "Centroid estimator (CONTRAfold)", marker = "p", linestyle = ":")
  line_3, = pyplot.plot(gammas, centroidfold_mccs, label = "CentroidFold (CONTRAfold)", marker = "s", linestyle = "--")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("MCC")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_ss_estimation.eps", bbox_inches = "tight")

def get_pos_neg_counts(params):
  (estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens) = params
  tp = tn = fp = fn = 0
  estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(estimated_ss_file_path))
  for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
    estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
    ref_ss, ref_flat_ss = ref_ss_and_flat_ss
    for i in range(0, rna_seq_len):
      for j in range(i + 1, rna_seq_len):
        estimated_bin = (i, j) in estimated_ss
        ref_bin = (i, j) in ref_ss
        if estimated_bin == ref_bin:
          if estimated_bin == True:
            tp += 1
          else:
            tn += 1
        else:
          if estimated_bin == True:
            fp += 1
          else:
            fn += 1
  return tp, tn, fp, fn

def final_sum(results):
  final_tp = final_tn = final_fp = final_fn = 0.
  for tp, tn, fp, fn in results:
    final_tp += tp
    final_tn += tn
    final_fp += fp
    final_fn += fn
  return (final_tp, final_tn, final_fp, final_fn)

def get_f1_score(ppv, sens):
  return 2 * ppv * sens / (ppv + sens)

def get_mcc(tp, tn, fp, fn):
  return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

def get_ppv(tp, fp):
  return tp / (tp + fp)

def get_sens(tp, fn):
  return tp / (tp + fn)

def get_fpr(tn, fp):
  return fp / (tn + fp)

if __name__ == "__main__":
  main()
