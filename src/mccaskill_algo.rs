use utils::*;

pub struct LogSsPpfMats {
  pub log_ss_ppf_mat: LogPpfMat,
  log_ss_ppf_mat_4_rightmost_base_pairings: LogPpfMat,
  log_ss_ppf_mat_4_base_pairings: LogPpfMat,
  log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls: LogPpfMat,
  log_ss_ppf_mat_4_rightmost_base_pairings_on_mls: LogPpfMat,
}

impl LogSsPpfMats {
  fn new(seq_len: usize) -> LogSsPpfMats {
    let ni_mat = vec![vec![NEG_INFINITY; seq_len]; seq_len];
    let mut log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls = ni_mat.clone();
    for i in 0 .. seq_len - 1 {
      log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[i + 1][i] = 0.;
    }
    LogSsPpfMats {
      log_ss_ppf_mat: vec![vec![0.; seq_len]; seq_len],
      log_ss_ppf_mat_4_rightmost_base_pairings: ni_mat.clone(),
      log_ss_ppf_mat_4_base_pairings: ni_mat.clone(),
      log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls: log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls,
      log_ss_ppf_mat_4_rightmost_base_pairings_on_mls: ni_mat,
    }
  }
}

pub const COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE: Energy = 0.;

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
        if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        let accessible_pp = (i + 1, j - 1);
        let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        for k in i + 2 .. j - 1 {
          let accessible_pp = (k, j - 1);
          let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
          if log_ss_ppf_4_base_pairing.is_finite() {
            let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
        for k in i + 2 .. j - 1 {
          let accessible_pp = (i + 1, k);
          let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
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
            let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
            if log_ss_ppf_4_base_pairing.is_finite() {
              let ep_of_term_4_log_pf = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
          }
        }
        if pp_closing_loop.1 - pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_ML {
          let mut eps_of_terms_4_log_coefficient = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_coefficient = 0.;
          eps_of_terms_4_log_coefficient.push(max_ep_of_term_4_log_coefficient);
          let tm = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
          let ml_tm_delta_fe = ML_TM_DELTA_FES[&(bp_closing_loop, tm)];
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * ml_tm_delta_fe;
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * FIVE_PRIME_DE_DELTA_FES[&(bp_closing_loop, seq[j - 1])];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * THREE_PRIME_DE_DELTA_FES[&(bp_closing_loop, seq[i + 1])];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
          let log_coefficient = logsumexp(&eps_of_terms_4_log_coefficient[..], max_ep_of_term_4_log_coefficient);
          for k in i + 1 .. j {
            let ep_of_term_4_log_pf = log_coefficient + log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[i + 1][k - 1] + log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings_on_mls[k][j - 1] - INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + if bp_closing_loop == AU || bp_closing_loop == UA || bp_closing_loop == GU || bp_closing_loop == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            if ep_of_term_4_log_pf.is_finite() {eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);}
          }
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in i + MIN_HL_LEN + 1 .. j + 1 {
        let accessible_pp = (i, k);
        let accessible_bp = invert_bp(&(seq[i], seq[k]));
        let log_ss_pf_4_bp = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
        if !log_ss_pf_4_bp.is_finite() {continue;}
        let log_coefficient = log_ss_pf_4_bp - INVERSE_TEMPERATURE * (if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.} + HELIX_INTERMOLECULAR_INIT_DELTA_FE);
        let ep_of_term_4_log_pf = log_coefficient;
        if ep_of_term_4_log_pf.is_finite() {
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        if i > 0 && k < seq_len - 1 {
          let ep_of_term_4_log_pf = log_coefficient - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(accessible_bp, invert_bp(&(seq[i - 1], seq[k + 1])))];
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
        if i > 0 {
          let ep_of_term_4_log_pf = log_coefficient - INVERSE_TEMPERATURE * FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])];
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
        if k < seq_len - 1 {
          let ep_of_term_4_log_pf = log_coefficient - INVERSE_TEMPERATURE * THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[k + 1])];
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = 0.;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      for k in i .. if j < MIN_HL_LEN {0} else {j - MIN_HL_LEN} {
        let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings[k][j];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let ep_of_term_4_log_pf = if k < 1 {0.} else {log_ss_ppf_mats.log_ss_ppf_mat[i][k - 1]} + log_ss_ppf_4_base_pairing;
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
      }
      log_ss_ppf_mats.log_ss_ppf_mat[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in i + MIN_HL_LEN + 1 .. j + 1 {
        let accessible_pp = (i, k);
        let accessible_bp = invert_bp(&(seq[accessible_pp.0], seq[accessible_pp.1]));
        let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
        if !log_ss_ppf_4_base_pairing.is_finite() {continue;}
        let log_coefficient = log_ss_ppf_4_base_pairing - INVERSE_TEMPERATURE * (COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (j - k) as Energy + if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.} + HELIX_INTERMOLECULAR_INIT_DELTA_FE);
        let ep_of_term_4_log_pf = log_coefficient;
        if ep_of_term_4_log_pf.is_finite() {
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        if i > 0 && k < seq_len - 1 {
          let ep_of_term_4_log_pf = log_coefficient - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(accessible_bp, invert_bp(&(seq[i - 1], seq[k + 1])))];
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
        if i > 0 {
          let ep_of_term_4_log_pf = log_coefficient - INVERSE_TEMPERATURE * FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])];
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
        if k < seq_len - 1 {
          let ep_of_term_4_log_pf = log_coefficient - INVERSE_TEMPERATURE * THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[k + 1])];
          if ep_of_term_4_log_pf.is_finite() {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings_on_mls[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in i .. if j < MIN_HL_LEN {0} else {j - MIN_HL_LEN} {
        let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
        let mut max_ep_of_sub_term_4_log_pf = -INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (k - i) as Energy;
        eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
        let ep_of_sub_term_4_log_pf = if k < 1 {0.} else {log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[i][k - 1]};
        if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
        eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
        let ep_of_term_4_log_pf = logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf) + log_ss_ppf_mats.log_ss_ppf_mat_4_rightmost_base_pairings_on_mls[k][j] - INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
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
      let accessible_pp = (i, j);
      let log_ss_ppf_4_base_pairing_1 = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[accessible_pp.0][accessible_pp.1];
      if log_ss_ppf_4_base_pairing_1.is_finite() {
        let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
        let mut eps_of_terms_4_log_coefficient = EpsOfTerms4LogProb::new();
        let mut max_ep_of_term_4_log_coefficient = 0. as LogProb;
        eps_of_terms_4_log_coefficient.push(max_ep_of_term_4_log_coefficient);
        let accessible_bp = invert_bp(&(seq[accessible_pp.0], seq[accessible_pp.1]));
        if i > 0 && j < seq_len - 1 {
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(accessible_bp, invert_bp(&(seq[i - 1], seq[j + 1])))];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
        }
        if i > 0 {
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
        }
        if j < seq_len - 1 {
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[j + 1])];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
        }
        let mut max_ep_of_term_4_log_prob = if accessible_pp.0 < 1 {0.} else {log_ss_ppf_mats.log_ss_ppf_mat[0][accessible_pp.0 - 1]} + log_ss_ppf_4_base_pairing_1 + if accessible_pp.1 > seq_len - 2 {0.} else {log_ss_ppf_mats.log_ss_ppf_mat[accessible_pp.1 + 1][seq_len - 1]} - log_ss_ppf + logsumexp(&eps_of_terms_4_log_coefficient[..], max_ep_of_term_4_log_coefficient) - INVERSE_TEMPERATURE * (if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.} + HELIX_INTERMOLECULAR_INIT_DELTA_FE);
        eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
        if i > 0 && j < seq_len - 1 {
          let pp_closing_loop = (i - 1, j + 1);
          let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
          if log_ss_ppf_4_base_pairing_2.is_finite() {
            let ep_of_term_4_log_prob = log_bpp_mat[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
        if i > 1 && j < seq_len - 1 {
          for k in 0 .. i - 1 {
            let pp_closing_loop = (k, j + 1);
            let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
            if log_ss_ppf_4_base_pairing_2.is_finite() {
              let ep_of_term_4_log_prob = log_bpp_mat[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
              if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            }
          }
        }
        if i > 0 && j < seq_len - 2 {
          for k in j + 2 .. seq_len {
            let pp_closing_loop = (i - 1, k);
            let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
            if log_ss_ppf_4_base_pairing_2.is_finite() {
              let ep_of_term_4_log_prob = log_bpp_mat[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
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
              let log_ss_ppf_4_base_pairing_2 = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
              if log_ss_ppf_4_base_pairing_2.is_finite() {
                let ep_of_term_4_log_prob = log_bpp_mat[pp_closing_loop.0][pp_closing_loop.1] + log_ss_ppf_4_base_pairing_1 - log_ss_ppf_4_base_pairing_2 - INVERSE_TEMPERATURE * get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as Energy;
                if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
              }
            }
          }
        }
        if i > 0 && j < seq_len - 1 {
          let accessible_bp = invert_bp(&(seq[accessible_pp.0], seq[accessible_pp.1]));
          let mut eps_of_terms_4_log_coefficient = EpsOfTerms4LogProb::new();
          let mut max_ep_of_term_4_log_coefficient = 0. as LogProb;
          eps_of_terms_4_log_coefficient.push(max_ep_of_term_4_log_coefficient);
          if i > 0 && j < seq_len - 1 {
            let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(accessible_bp, invert_bp(&(seq[i - 1], seq[j + 1])))];
            if ep_of_term_4_log_coefficient.is_finite() {
              if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
              eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
            }
          }
          if i > 0 {
            let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])];
            if ep_of_term_4_log_coefficient.is_finite() {
              if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
              eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
            }
          }
          if j < seq_len - 1 {
            let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[j + 1])];
            if ep_of_term_4_log_coefficient.is_finite() {
              if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
              eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
            }
          }
          let log_coefficient = log_ss_ppf_4_base_pairing_1 - INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.} + HELIX_INTERMOLECULAR_INIT_DELTA_FE) + logsumexp(&eps_of_terms_4_log_coefficient[..], max_ep_of_term_4_log_coefficient);
          for k in 0 .. i {
            let mut eps_of_sub_terms_4_log_prob = EpsOfTerms4LogProb::new();
            let log_ss_ppf_4_at_least_1_base_pairings_on_mls = log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[k + 1][i - 1];
            let mut max_ep_of_sub_term_4_log_prob = log_prob_mat_4_mls_2[k][j] + log_ss_ppf_4_at_least_1_base_pairings_on_mls;
            eps_of_sub_terms_4_log_prob.push(max_ep_of_sub_term_4_log_prob);
            let log_prob_4_mls = log_prob_mat_4_mls_1[k][j];
            let ep_of_sub_term_4_log_prob = log_prob_4_mls - INVERSE_TEMPERATURE * COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (i - k - 1) as Energy;
            if max_ep_of_sub_term_4_log_prob < ep_of_sub_term_4_log_prob {max_ep_of_sub_term_4_log_prob = ep_of_sub_term_4_log_prob;}
            eps_of_sub_terms_4_log_prob.push(ep_of_sub_term_4_log_prob);
            let ep_of_sub_term_4_log_prob = log_prob_4_mls + log_ss_ppf_4_at_least_1_base_pairings_on_mls;
            if max_ep_of_sub_term_4_log_prob < ep_of_sub_term_4_log_prob {max_ep_of_sub_term_4_log_prob = ep_of_sub_term_4_log_prob;}
            eps_of_sub_terms_4_log_prob.push(ep_of_sub_term_4_log_prob);
            let ep_of_term_4_log_prob = log_coefficient + logsumexp(&eps_of_sub_terms_4_log_prob[..], max_ep_of_sub_term_4_log_prob);
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
        if eps_of_terms_4_log_prob.len() > 0 {
          log_bpp_mat[accessible_pp.0][accessible_pp.1] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
        }
      }
      let mut eps_of_terms_4_log_prob_1 = EpsOfTerms4LogProb::new();
      let mut eps_of_terms_4_log_prob_2 = EpsOfTerms4LogProb::new();
      let mut max_ep_of_term_4_log_prob_1 = NEG_INFINITY;
      let mut max_ep_of_term_4_log_prob_2 = NEG_INFINITY;
      for k in j + 1 .. seq_len {
        let pp_closing_loop = (i, k);
        let log_ss_ppf_4_base_pairing = log_ss_ppf_mats.log_ss_ppf_mat_4_base_pairings[pp_closing_loop.0][pp_closing_loop.1];
        if log_ss_ppf_4_base_pairing.is_finite() {
          let bp_closing_loop = (seq[i], seq[k]);
          let mut eps_of_terms_4_log_coefficient = EpsOfTerms4LogProb::new();
          let mut max_ep_of_term_4_log_coefficient = 0. as LogProb;
          eps_of_terms_4_log_coefficient.push(max_ep_of_term_4_log_coefficient);
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(bp_closing_loop, (seq[i + 1], seq[j - 1]))];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * FIVE_PRIME_DE_DELTA_FES[&(bp_closing_loop, seq[k - 1])];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
          let ep_of_term_4_log_coefficient = -INVERSE_TEMPERATURE * THREE_PRIME_DE_DELTA_FES[&(bp_closing_loop, seq[i + 1])];
          if ep_of_term_4_log_coefficient.is_finite() {
            if max_ep_of_term_4_log_coefficient < ep_of_term_4_log_coefficient {max_ep_of_term_4_log_coefficient = ep_of_term_4_log_coefficient;}
            eps_of_terms_4_log_coefficient.push(ep_of_term_4_log_coefficient);
          }
          let log_bpp = log_bpp_mat[pp_closing_loop.0][pp_closing_loop.1];
          let helix_au_or_gu_end_penalty_delta_fe = if bp_closing_loop == AU || bp_closing_loop == UA || bp_closing_loop == GU || bp_closing_loop == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
          let ep_of_term_4_log_prob = log_bpp - log_ss_ppf_4_base_pairing + log_ss_ppf_mats.log_ss_ppf_mat_4_at_least_1_base_pairings_on_mls[j + 1][pp_closing_loop.1 - 1] + logsumexp(&eps_of_terms_4_log_coefficient[..], max_ep_of_term_4_log_coefficient) - INVERSE_TEMPERATURE * helix_au_or_gu_end_penalty_delta_fe;
          if max_ep_of_term_4_log_prob_1 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_1 = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob_1.push(ep_of_term_4_log_prob);
          let ep_of_term_4_log_prob = log_bpp - log_ss_ppf_4_base_pairing + logsumexp(&eps_of_terms_4_log_coefficient[..], max_ep_of_term_4_log_coefficient) - INVERSE_TEMPERATURE * (COEFFICENT_4_TERM_OF_NUM_OF_UNPAIRED_BASES_ON_INIT_ML_DELTA_FE * (pp_closing_loop.1 - j - 1) as Energy + helix_au_or_gu_end_penalty_delta_fe);
          if max_ep_of_term_4_log_prob_2 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_2 = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob_2.push(ep_of_term_4_log_prob);
        }
      }
      if eps_of_terms_4_log_prob_1.len() > 0 {
        log_prob_mat_4_mls_1[i][j] = logsumexp(&eps_of_terms_4_log_prob_1[..], max_ep_of_term_4_log_prob_1);
        log_prob_mat_4_mls_2[i][j] = logsumexp(&eps_of_terms_4_log_prob_2[..], max_ep_of_term_4_log_prob_2);
      }
    }
  }
  log_bpp_mat
}
