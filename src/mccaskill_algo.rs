use utils::*;

pub struct LogSsPpfMatrices {
  pub log_ss_ppf_matrix: LogPpfMatrix,
  log_ss_ppf_matrix_4_rightmost_base_pairings: LogPpfMatrix,
  log_ss_ppf_matrix_4_base_pairings: LogPpfMatrix,
  log_ss_ppf_matrix_4_at_least_1_base_pairings_on_mls: LogPpfMatrix,
  log_ss_ppf_matrix_4_rightmost_base_pairings_on_mls: LogPpfMatrix,
}

impl LogSsPpfMatrices {
  fn new(seq_len: usize) -> LogSsPpfMatrices {
    let ni_matrix = vec![vec![NEG_INFINITY; seq_len]; seq_len];
    LogSsPpfMatrices {
      log_ss_ppf_matrix: vec![vec![0.; seq_len]; seq_len],
      log_ss_ppf_matrix_4_rightmost_base_pairings: ni_matrix.clone(),
      log_ss_ppf_matrix_4_base_pairings: ni_matrix.clone(),
      log_ss_ppf_matrix_4_at_least_1_base_pairings_on_mls: ni_matrix.clone(),
      log_ss_ppf_matrix_4_rightmost_base_pairings_on_mls: ni_matrix,
    }
  }
}

pub const CONST_4_INIT_ML_DELTA_FE: Energy = 10.1;
pub const COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE: Energy = -0.3;
pub const COEFFICENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE: Energy = -0.3;

#[inline]
pub fn mccaskill_algo(seq: SeqSlice) -> ProbMatrix {
  let seq_len = seq.len();
  let log_ss_ppf_matrices = get_log_ss_ppf_matrices(&seq[..], seq_len);
  let log_bpp_matrix = get_log_base_pairing_prob_matrix(&seq[..], &log_ss_ppf_matrices, seq_len);
  get_bpp_matrix(&log_bpp_matrix)
}

#[inline]
fn get_bpp_matrix(log_bpp_matrix: &LogProbMatrix) -> ProbMatrix {
  log_bpp_matrix.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect()
}

#[inline]
pub fn get_bpp_matrix_and_unpaired_base_probs(seq: SeqSlice) -> (ProbMatrix, Probs) {
  let seq_len = seq.len();
  let log_ss_ppf_matrices = get_log_ss_ppf_matrices(&seq[..], seq_len);
  let log_bpp_matrix = get_log_base_pairing_prob_matrix(&seq[..], &log_ss_ppf_matrices, seq_len);
  let mut ubps = vec![NEG_INFINITY; seq_len];
  for i in 0 .. seq_len {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = NEG_INFINITY;
    for j in 0 .. i {
      if j == i {continue;}
      let pp = if j < i {(j, i)} else {(i, j)};
      let ep_of_term_4_log_prob = log_bpp_matrix[pp.0][pp.1];
      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
    }
    ubps[i] = 1. - logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob).exp();
  }
  (get_bpp_matrix(&log_bpp_matrix), ubps)
}

#[inline]
pub fn get_log_ss_ppf_matrices(seq: SeqSlice, seq_len: usize) -> LogSsPpfMatrices {
  let mut log_ss_ppf_matrices = LogSsPpfMatrices::new(seq_len);
  for sub_seq_len in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. seq_len + 1 {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let pp_closing_loop = (i, j);
      let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if pp_closing_loop.1 - pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL && is_canonical(&bp_closing_loop) {
        let ep_of_term_4_log_pf = -INVERSE_TEMPERATURE * get_hl_fe(seq, &pp_closing_loop) as Energy;
        if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        let accessible_pp = (i + 1, j - 1);
        let log_ss_ppf_4_base_pairing = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[accessible_pp.0][accessible_pp.1];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        for k in i + 2 .. j - 1 {
          let accessible_pp = (k, j - 1);
          let log_ss_ppf_4_base_pairing = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[accessible_pp.0][accessible_pp.1];
          if log_ss_ppf_4_base_pairing.is_finite() {
            let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
        for k in i + 2 .. j - 1 {
          let accessible_pp = (i + 1, k);
          let log_ss_ppf_4_base_pairing = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[accessible_pp.0][accessible_pp.1];
          if log_ss_ppf_4_base_pairing.is_finite() {
            let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
        for k in i + 2 .. min(i + MAX_IL_LEN, j - 2) {
          let min_l = j as isize - (MAX_IL_LEN as isize - (k - i - 1) as isize) - 1;
          let min_l = if min_l > k as isize + 1 {min_l as usize} else {k + 1};
          for l in min_l .. j - 1 {
            debug_assert!(1 < j - l - 1 + k - i - 1 && j - l - 1 + k - i - 1 <= MAX_IL_LEN);
            let accessible_pp = (k, l);
            let log_ss_ppf_4_base_pairing = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[accessible_pp.0][accessible_pp.1];
            if log_ss_ppf_4_base_pairing.is_finite() {
              let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
          }
        }
        if pp_closing_loop.1 - pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_ML {
          for k in i + 1 .. j {
            let ep_of_term_4_log_pf = log_ss_ppf_matrices.log_ss_ppf_matrix_4_at_least_1_base_pairings_on_mls[i + 1][k - 1] + log_ss_ppf_matrices.log_ss_ppf_matrix_4_rightmost_base_pairings_on_mls[k][j - 1] - INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            if ep_of_term_4_log_pf.is_finite() {eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);}
          }
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in i + MIN_HL_LEN + 1 .. j + 1 {
        let accessible_pp = (i, k);
        let ep_of_term_4_log_pf = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[accessible_pp.0][accessible_pp.1];
        if ep_of_term_4_log_pf.is_finite() {
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_matrices.log_ss_ppf_matrix_4_rightmost_base_pairings[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = 0.;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      for k in i .. if j < MIN_HL_LEN {0} else {j - MIN_HL_LEN} {
        let log_ss_ppf_4_base_pairing = log_ss_ppf_matrices.log_ss_ppf_matrix_4_rightmost_base_pairings[k][j];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let ep_of_term_4_log_pf = if k < 1 {0.} else {log_ss_ppf_matrices.log_ss_ppf_matrix[i][k - 1]} + log_ss_ppf_4_base_pairing;
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_matrices.log_ss_ppf_matrix[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      let min_k = i + 1;
      let sup_k = j + 1;
      for k in min_k .. sup_k {
        let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
        let mut max_ep_of_sub_term_4_log_pf = -INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (k - i) as Energy;
        eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
        let ep_of_sub_term_4_log_pf = log_ss_ppf_matrices.log_ss_ppf_matrix_4_at_least_1_base_pairings_on_mls[i][k - 1];
        if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
        eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
        let ep_of_term_4_log_pf = logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf) + log_ss_ppf_matrices.log_ss_ppf_matrix_4_rightmost_base_pairings_on_mls[k][j] - INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
        if ep_of_term_4_log_pf.is_finite() {eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);}
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_matrices.log_ss_ppf_matrix_4_at_least_1_base_pairings_on_mls[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in min_k .. sup_k {
        let accessible_pp = (i, k);
        let log_ss_ppf_4_base_pairing = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[accessible_pp.0][accessible_pp.1];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (j - k) as Energy;
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_matrices.log_ss_ppf_matrix_4_rightmost_base_pairings_on_mls[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
    }
  }
  log_ss_ppf_matrices
}

#[inline]
fn get_log_base_pairing_prob_matrix(seq: SeqSlice, log_ss_ppf_matrices: &LogSsPpfMatrices, seq_len: usize) -> LogProbMatrix {
  let log_ss_ppf = log_ss_ppf_matrices.log_ss_ppf_matrix[0][seq_len - 1];
  let mut log_bpp_matrix = vec![vec![NEG_INFINITY; seq_len]; seq_len];
  let mut log_prob_matrix_4_mls_1 = log_bpp_matrix.clone();
  let mut log_prob_matrix_4_mls_2 = log_bpp_matrix.clone();
  for sub_seq_len in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. seq_len + 1).rev() {
    for i in 0 .. seq_len - sub_seq_len {
      let j = i + sub_seq_len - 1;
      let accessible_pp = (i, j);
      let log_ss_ppf_4_base_pairing_1 = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[accessible_pp.0][accessible_pp.1];
      if log_ss_ppf_4_base_pairing_1.is_finite() {
        let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
        let mut max_ep_of_term_4_log_prob = if accessible_pp.0 < 1 {0.} else {log_ss_ppf_matrices.log_ss_ppf_matrix[0][accessible_pp.0 - 1]} + log_ss_ppf_4_base_pairing_1 + if accessible_pp.1 > seq_len - 2 {0.} else {log_ss_ppf_matrices.log_ss_ppf_matrix[accessible_pp.1 + 1][seq_len - 1]} - log_ss_ppf;
        eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
        if i > 0 && j < seq_len - 1 {
          let pp_closing_loop = (i - 1, j + 1);
          let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
          if log_ss_ppf_4_base_pairing_2.is_finite() {
            let ep_of_term_4_log_prob = log_bpp_matrix[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
        if i > 1 && j < seq_len - 1 {
          for k in 0 .. i - 1 {
            let pp_closing_loop = (k, j + 1);
            let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
            if log_ss_ppf_4_base_pairing_2.is_finite() {
              let ep_of_term_4_log_prob = log_bpp_matrix[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
              if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            }
          }
        }
        if i > 0 && j < seq_len - 2 {
          for k in j + 2 .. seq_len {
            let pp_closing_loop = (i - 1, k);
            let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
            if log_ss_ppf_4_base_pairing_2.is_finite() {
              let ep_of_term_4_log_prob = log_bpp_matrix[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
              if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            }
          }
        }
        if i > 1 && j < seq_len - 2 {
          let min_k = i as isize - MAX_IL_LEN as isize;
          let min_k = if min_k > 1 {min_k as usize} else {1};
          for k in min_k .. i - 1 {
            let sup_l = j as isize + (MAX_IL_LEN as isize - (i - k - 1) as isize) + 2;
            let sup_l = if sup_l >= 0 && sup_l < seq_len as isize {sup_l as usize} else if sup_l < 0 {1} else {seq_len};
            for l in j + 2 .. sup_l {
              debug_assert!(1 < l - j - 1 + i - k - 1 && l - j - 1 + i - k - 1 <= MAX_IL_LEN);
              let pp_closing_loop = (k, l);
              let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
              if log_ss_ppf_4_base_pairing_2.is_finite() {
                let ep_of_term_4_log_prob = log_bpp_matrix[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
                if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
              }
            }
          }
        }
        if i > 0 && j < seq_len - 1 {
          for k in 0 .. i {
            let mut eps_of_sub_terms_4_log_prob = EpsOfTerms4LogProb::new();
            let log_ss_ppf_4_at_least_1_base_pairings_on_mls = log_ss_ppf_matrices.log_ss_ppf_matrix_4_at_least_1_base_pairings_on_mls[k + 1][i - 1];
            let mut max_ep_of_sub_term_4_log_prob = log_prob_matrix_4_mls_2[k][j] + log_ss_ppf_4_at_least_1_base_pairings_on_mls;
            eps_of_sub_terms_4_log_prob.push(max_ep_of_sub_term_4_log_prob);
            let log_prob_4_mls = log_prob_matrix_4_mls_1[k][j];
            let ep_of_sub_term_4_log_prob = log_prob_4_mls - INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (i - k - 1) as Energy;
            if max_ep_of_sub_term_4_log_prob < ep_of_sub_term_4_log_prob {max_ep_of_sub_term_4_log_prob = ep_of_sub_term_4_log_prob;}
            eps_of_sub_terms_4_log_prob.push(ep_of_sub_term_4_log_prob);
            let ep_of_sub_term_4_log_prob = log_prob_4_mls + log_ss_ppf_4_at_least_1_base_pairings_on_mls;
            if max_ep_of_sub_term_4_log_prob < ep_of_sub_term_4_log_prob {max_ep_of_sub_term_4_log_prob = ep_of_sub_term_4_log_prob;}
            eps_of_sub_terms_4_log_prob.push(ep_of_sub_term_4_log_prob);
            let ep_of_term_4_log_prob = log_ss_ppf_4_base_pairing_1 - INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) + logsumexp(&eps_of_sub_terms_4_log_prob[..], max_ep_of_sub_term_4_log_prob);
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
        if eps_of_terms_4_log_prob.len() > 0 {
          log_bpp_matrix[accessible_pp.0][accessible_pp.1] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
        }
      }
      let mut eps_of_terms_4_log_prob_1 = EpsOfTerms4LogProb::new();
      let mut eps_of_terms_4_log_prob_2 = EpsOfTerms4LogProb::new();
      let mut max_ep_of_term_4_log_prob_1 = NEG_INFINITY;
      let mut max_ep_of_term_4_log_prob_2 = NEG_INFINITY;
      for k in j + 1 .. seq_len {
        let pp_closing_loop = (i, k);
        let log_ss_ppf_4_base_pairing = log_ss_ppf_matrices.log_ss_ppf_matrix_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let log_bpp = log_bpp_matrix[pp_closing_loop.0][pp_closing_loop.1];
          let ep_of_term_4_log_prob = log_bpp - log_ss_ppf_4_base_pairing + log_ss_ppf_matrices.log_ss_ppf_matrix_4_at_least_1_base_pairings_on_mls[j + 1][pp_closing_loop.1 - 1];
          if max_ep_of_term_4_log_prob_1 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_1 = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob_1.push(ep_of_term_4_log_prob);
          let ep_of_term_4_log_prob = log_bpp - log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (pp_closing_loop.1 - j - 1) as Energy;
          if max_ep_of_term_4_log_prob_2 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_2 = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob_2.push(ep_of_term_4_log_prob);
        }
      }
      if eps_of_terms_4_log_prob_1.len() > 0 {
        log_prob_matrix_4_mls_1[i][j] = logsumexp(&eps_of_terms_4_log_prob_1[..], max_ep_of_term_4_log_prob_1);
        log_prob_matrix_4_mls_2[i][j] = logsumexp(&eps_of_terms_4_log_prob_2[..], max_ep_of_term_4_log_prob_2);
      }
    }
  }
  log_bpp_matrix
}
