use utils::*;

pub struct LogSsPpfMats {
  pub log_ss_ppf_mat: LogPpfMat,
  log_ss_ppf_mat_4_rightmost_base_pairings: LogPpfMat,
  log_ss_ppf_mat_4_base_pairings: LogPpfMat,
  log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls: LogPpfMat,
}

impl LogSsPpfMats {
  fn new(seq_len: usize) -> LogSsPpfMats {
    let ni_mat = vec![vec![NEG_INFINITY; seq_len]; seq_len];
    LogSsPpfMats {
      log_ss_ppf_mat: vec![vec![0.; seq_len]; seq_len],
      log_ss_ppf_mat_4_rightmost_base_pairings: ni_mat.clone(),
      log_ss_ppf_mat_4_base_pairings: ni_mat.clone(),
      log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls: ni_mat,
    }
  }
}

pub const CONST_4_INIT_ML_DELTA_FE: FreeEnergy = 9.3;
pub const COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE: FreeEnergy = -0.9;

#[inline]
pub fn mccaskill_algo(seq: SeqSlice) -> ProbMat {
  let seq_len = seq.len();
  let log_ss_ppf_mats = get_log_ss_ppf_mats(seq, seq_len);
  let log_bpp_mat = get_log_base_pairing_prob_mat(seq, &log_ss_ppf_mats, seq_len);
  get_bpp_mat(&log_bpp_mat)
}

#[inline]
pub fn get_bpp_mat(log_bpp_mat: &LogProbMat) -> ProbMat {
  log_bpp_mat.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect()
}

#[inline]
pub fn get_log_bpp_mat(seq: SeqSlice) -> LogProbMat {
  let seq_len = seq.len();
  let log_ss_ppf_mats = get_log_ss_ppf_mats(&seq[..], seq_len);
  get_log_base_pairing_prob_mat(&seq[..], &log_ss_ppf_mats, seq_len)
}

#[inline]
pub fn get_bpp_mat_and_nbpps(seq: SeqSlice) -> (ProbMat, Probs) {
  let seq_len = seq.len();
  let log_ss_ppf_mats = get_log_ss_ppf_mats(&seq[..], seq_len);
  let log_bpp_mat = get_log_base_pairing_prob_mat(&seq[..], &log_ss_ppf_mats, seq_len);
  let mut nbpps = vec![NEG_INFINITY; seq_len];
  for i in 0 .. seq_len {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = NEG_INFINITY;
    for j in 0 .. seq_len {
      if j == i {continue;}
      let pp = if j < i {(j, i)} else {(i, j)};
      let ep_of_term_4_log_prob = log_bpp_mat[pp.0][pp.1];
      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
    }
    nbpps[i] = 1. - logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob).exp();
  }
  (get_bpp_mat(&log_bpp_mat), nbpps)
}

#[inline]
pub fn get_log_ss_ppf_mats(seq: SeqSlice, seq_len: usize) -> LogSsPpfMats {
  let mut log_ss_ppf_mats = LogSsPpfMats::new(seq_len);
  for sub_seq_len in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. seq_len + 1 {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let pp_closing_loop = (i, j);
      let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if pp_closing_loop.1 - pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL && is_canonical(&bp_closing_loop) {
        let ep_of_term_4_log_pf = -INVERSE_TEMPERATURE * get_hl_fe(seq, &pp_closing_loop) as Energy;
        if ep_of_term_4_log_pf.is_finite() {
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        for k in i + 1 .. j - 1 {
          for l in k + 1 .. j {
            if j - l - 1 + k - i - 1 > MAX_2_LOOP_LEN {continue;}
            let accessible_pp = (k, l);
            let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
            if log_ss_ppf_4_base_pairing.is_finite() {
              let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
          }
        }
        for k in i + 1 .. j {
          let ep_of_term_4_log_pf = log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[i + 1][k - 1] + log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings[k][j - 1] - INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop), invert_bp(&(seq[i + 1], seq[j - 1])))] + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in i + 1 .. j + 1 {
        let accessible_pp = (i, k);
        let accessible_bp = (seq[i], seq[k]);
        let log_ss_pf_4_bp = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
        if !log_ss_pf_4_bp.is_finite() {continue;}
        let ep_of_term_4_log_pf = log_ss_pf_4_bp - INVERSE_TEMPERATURE * (if i > 0 && k < seq_len - 1 {
          ML_TM_DELTA_FES[&(accessible_bp, (seq[i - 1], seq[k + 1]))]
        } else if i > 0 {
          FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])]
        } else if k < seq_len - 1 {
          THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[k + 1])]
        } else {
          0.
        } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
        if ep_of_term_4_log_pf.is_finite() {
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = 0.;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      for k in i .. j {
        let log_ss_ppf_4_rightmost_base_pairings = log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings[k][j];
        if log_ss_ppf_4_rightmost_base_pairings.is_finite() {
          let ep_of_term_4_log_pf = if i == 0 && k == 0 {0.} else {log_ss_ppf_mats.log_ss_ppf_mat[i][k - 1]} + log_ss_ppf_4_rightmost_base_pairings;
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
      }
      if eps_of_terms_4_log_pf.len() > 1 {
        log_ss_ppf_mats.log_ss_ppf_mat[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings[i][j] - INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
      if max_ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      }
      for k in i + 1 .. j {
        let log_ss_ppf_4_rightmost_base_pairings = log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings[k][j];
        let ep_of_term_4_log_pf = log_ss_ppf_4_rightmost_base_pairings - INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        if ep_of_term_4_log_pf.is_finite() {
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        let ep_of_term_4_log_pf = log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[i][k - 1] + log_ss_ppf_4_rightmost_base_pairings - INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        if ep_of_term_4_log_pf.is_finite() {
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
    }
  }
  log_ss_ppf_mats
}

#[inline]
fn get_log_base_pairing_prob_mat(seq: SeqSlice, log_ss_ppf_mats: &LogSsPpfMats, seq_len: usize) -> LogProbMat {
  let log_ss_ppf = log_ss_ppf_mats.log_ss_ppf_mat[0][seq_len - 1];
  let mut log_bpp_mat = vec![vec![NEG_INFINITY; seq_len]; seq_len];
  let mut log_prob_mat_4_mls_1 = log_bpp_mat.clone();
  let mut log_prob_mat_4_mls_2 = log_bpp_mat.clone();
  for sub_seq_len in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. seq_len + 1).rev() {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let mut eps_of_terms_4_log_prob_1 = EpsOfTerms4LogProb::new();
      let mut eps_of_terms_4_log_prob_2 = EpsOfTerms4LogProb::new();
      let mut max_ep_of_term_4_log_prob_1 = NEG_INFINITY;
      let mut max_ep_of_term_4_log_prob_2 = NEG_INFINITY;
      for k in j + 1 .. seq_len {
        let pp_closing_loop = (i, k);
        let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let log_bpp = log_bpp_mat[pp_closing_loop.0][pp_closing_loop.1];
          let bp_closing_loop = (seq[i], seq[k]);
          let ml_tm_delta_fe = ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop), invert_bp(&(seq[i + 1], seq[k - 1])))];
          let log_coefficient = log_bpp - INVERSE_TEMPERATURE * (ml_tm_delta_fe + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
          let ep_of_term_4_log_prob = log_coefficient - log_ss_ppf_4_base_pairing + log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[j + 1][pp_closing_loop.1 - 1];
          if ep_of_term_4_log_prob.is_finite() {
            if max_ep_of_term_4_log_prob_1 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_1 = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob_1.push(ep_of_term_4_log_prob);
          }
          let ep_of_term_4_log_prob = log_coefficient - log_ss_ppf_4_base_pairing;
          if ep_of_term_4_log_prob.is_finite() {
            if max_ep_of_term_4_log_prob_2 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_2 = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob_2.push(ep_of_term_4_log_prob);
          }
        }
      }
      if eps_of_terms_4_log_prob_1.len() > 0 {
        log_prob_mat_4_mls_1[i][j] = logsumexp(&eps_of_terms_4_log_prob_1[..], max_ep_of_term_4_log_prob_1);
      }
      if eps_of_terms_4_log_prob_2.len() > 0 {
        log_prob_mat_4_mls_2[i][j] = logsumexp(&eps_of_terms_4_log_prob_2[..], max_ep_of_term_4_log_prob_2);
      }
      let accessible_pp = (i, j);
      let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
      let log_ss_ppf_4_base_pairing_1 = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
      if !log_ss_ppf_4_base_pairing_1.is_finite() {continue;}
      let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
      let mut max_ep_of_term_4_log_prob = if accessible_pp.0 < 1 {0.} else {log_ss_ppf_mats.log_ss_ppf_mat[0][accessible_pp.0 - 1]} + log_ss_ppf_4_base_pairing_1 + if accessible_pp.1 > seq_len - 2 {0.} else {log_ss_ppf_mats.log_ss_ppf_mat[accessible_pp.1 + 1][seq_len - 1]} - log_ss_ppf - INVERSE_TEMPERATURE * (if i > 0 && j < seq_len - 1 {
        ML_TM_DELTA_FES[&(accessible_bp, (seq[i - 1], seq[j + 1]))]
      } else if i > 0 {
        FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])]
      } else if j < seq_len - 1 {
        THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[j + 1])]
      } else {
        0.
      } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
      if max_ep_of_term_4_log_prob.is_finite() {
        eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
      }
      for k in 0 .. i {
        for l in j + 1 .. seq_len {
          if l - j - 1 + i - k - 1 > MAX_2_LOOP_LEN {continue;}
          let pp_closing_loop = (k, l);
          let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
          if log_ss_ppf_4_base_pairing_2.is_finite() {
            let ep_of_term_4_log_prob = log_bpp_mat[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
            if ep_of_term_4_log_prob.is_finite() {
              if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            }
          }
        }
      }
      let log_coefficient = log_ss_ppf_4_base_pairing_1 - INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + if i > 0 && j < seq_len - 1 {
        ML_TM_DELTA_FES[&(accessible_bp, (seq[i - 1], seq[j + 1]))]
      } else if i > 0 {
        FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])]
      } else if j < seq_len - 1 {
        THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[j + 1])]
      } else {
        0.
      } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
      for k in 0 .. i {
        let log_ss_ppf_4_at_least_1_base_pairings_on_mls = log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[k + 1][i - 1];
        let ep_of_term_4_log_prob = log_coefficient + log_prob_mat_4_mls_2[k][j] + log_ss_ppf_4_at_least_1_base_pairings_on_mls;
        if ep_of_term_4_log_prob.is_finite() {
          if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
        }
        let log_prob_4_mls = log_prob_mat_4_mls_1[k][j];
        let ep_of_term_4_log_prob = log_coefficient + log_prob_4_mls;
        if ep_of_term_4_log_prob.is_finite() {
          if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
        }
        let ep_of_term_4_log_prob = log_coefficient + log_prob_4_mls + log_ss_ppf_4_at_least_1_base_pairings_on_mls;
        if ep_of_term_4_log_prob.is_finite() {
          if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
        }
      }
      if eps_of_terms_4_log_prob.len() > 0 {
        log_bpp_mat[accessible_pp.0][accessible_pp.1] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
      }
    }
  }
  log_bpp_mat
}
