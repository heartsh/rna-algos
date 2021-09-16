#! /usr/bin/env python

import compile_rna_fams
import run_ss_estimation_programs
import get_stats_of_ss_estimation_programs

def main():
  compile_rna_fams.main()
  run_ss_estimation_programs.main()
  get_stats_of_ss_estimation_programs.main()

if __name__ == "__main__":
  main()
