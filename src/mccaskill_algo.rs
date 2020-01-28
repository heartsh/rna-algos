use utils::*;

pub struct SsPartFuncMats {
  pub part_func_flag_mat: PartFuncFlagMat,
  pub part_func_mat_4_rightmost_base_pairings: PartFuncMat,
  pub part_func_mat_4_base_pairings: SparsePartFuncMat,
  pub part_func_mat_4_at_least_1_base_pairings_on_mls: PartFuncMat,
}
pub struct SsMaxFreeEnergyMats {
  pub max_free_energy_mat: FreeEnergyMat,
  pub max_free_energy_mat_4_rightmost_base_pairings: FreeEnergyMat,
  pub max_free_energy_mat_4_base_pairings: SparseFreeEnergyMat,
  pub max_free_energy_mat_4_at_least_1_base_pairings_on_mls: FreeEnergyMat,
}
pub type FreeEnergies = Vec<FreeEnergy>;
pub type FreeEnergyMat = Vec<FreeEnergies>;
pub type SparseFreeEnergyMat = HashMap<PosPair, FreeEnergy, Hasher>;
#[derive(Clone)]
pub struct PartFuncFlag {
  pub part_func: PartFunc,
  pub flag: bool,
}
pub type PartFuncFlags = Vec<PartFuncFlag>;
pub type PartFuncFlagMat = Vec<PartFuncFlags>;
pub type SparsePartFuncMat = HashMap<PosPair, PartFunc, Hasher>;
pub type SparseProbMat = HashMap<PosPair, Prob, Hasher>;

impl SsPartFuncMats {
  fn new(seq_len: usize) -> SsPartFuncMats {
    let zero_mat = vec![vec![0.; seq_len]; seq_len];
    SsPartFuncMats {
      part_func_flag_mat: vec![vec![PartFuncFlag::new(); seq_len]; seq_len],
      part_func_mat_4_rightmost_base_pairings: zero_mat.clone(),
      part_func_mat_4_base_pairings: SparsePartFuncMat::default(),
      part_func_mat_4_at_least_1_base_pairings_on_mls: zero_mat,
    }
  }
}

impl PartFuncFlag {
  fn new() -> PartFuncFlag {
    PartFuncFlag {
      part_func: 1.,
      flag: false,
    }
  }
}

impl SsMaxFreeEnergyMats {
  fn new(seq_len: usize) -> SsMaxFreeEnergyMats {
    let ni_mat = vec![vec![NEG_INFINITY; seq_len]; seq_len];
    SsMaxFreeEnergyMats {
      max_free_energy_mat: vec![vec![0.; seq_len]; seq_len],
      max_free_energy_mat_4_rightmost_base_pairings: ni_mat.clone(),
      max_free_energy_mat_4_base_pairings: SparseFreeEnergyMat::default(),
      max_free_energy_mat_4_at_least_1_base_pairings_on_mls: ni_mat,
    }
  }
}

pub const CONST_4_INIT_ML_DELTA_FE: FreeEnergy = - INVERSE_TEMPERATURE * 9.3;
pub const COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE: FreeEnergy = - INVERSE_TEMPERATURE * (-0.9);
lazy_static! {
  pub static ref EXP_CONST_4_INIT_ML_DELTA_FE: FreeEnergy = CONST_4_INIT_ML_DELTA_FE.exp();
  pub static ref EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE: FreeEnergy = COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE.exp();
}

#[inline]
pub fn mccaskill_algo(seq: SeqSlice) -> (SparseProbMat, FreeEnergy) {
  let seq_len = seq.len();
  let max_free_energy = get_max_free_energy(seq, seq_len);
  let invert_exp_max_free_energy = 1. / max_free_energy.exp();
  let ss_part_func_mats = get_ss_part_func_mats(seq, seq_len, invert_exp_max_free_energy);
  let bpp_mat = get_base_pairing_prob_mat(seq, &ss_part_func_mats, seq_len, invert_exp_max_free_energy);
  (bpp_mat, max_free_energy)
}

#[inline]
pub fn get_bpp_and_unpair_prob_mats(seq: SeqSlice) -> (SparseProbMat, Probs, FreeEnergy) {
  let seq_len = seq.len();
  let (bpp_mat, max_free_energy) = mccaskill_algo(&seq[..]);
  let mut unpair_prob_mat = vec![1.; seq_len];
  for i in 0 .. seq_len {
    let mut sum = 0.;
    for j in 0 .. seq_len {
      if j == i {continue;}
      let pp = if j < i {(j, i)} else {(i, j)};
      if !bpp_mat.contains_key(&pp) {continue;}
      sum += bpp_mat[&pp];
    }
    unpair_prob_mat[i] = 1. - sum;
  }
  (bpp_mat, unpair_prob_mat, max_free_energy)
}

#[inline]
pub fn get_max_free_energy(seq: SeqSlice, seq_len: usize) -> FreeEnergy {
  let mut ss_max_free_energy_mats = SsMaxFreeEnergyMats::new(seq_len);
  for sub_seq_len in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. seq_len + 1 {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let pp_closing_loop = (i, j);
      let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
      let mut max_free_energy;
      if pp_closing_loop.1 - pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL && is_canonical(&bp_closing_loop) {
        max_free_energy = get_hl_fe(seq, &pp_closing_loop);
        for k in i + 1 .. j - 1 {
          for l in k + 1 .. j {
            if j - l - 1 + k - i - 1 > MAX_2_LOOP_LEN {continue;}
            let accessible_pp = (k, l);
            if !ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
            let ss_max_free_energy_4_base_pairing = ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings[&accessible_pp];
            let twoloop_free_energy = ss_max_free_energy_4_base_pairing + get_2_loop_fe(seq, &pp_closing_loop, &accessible_pp);
            if twoloop_free_energy > max_free_energy {max_free_energy = twoloop_free_energy};
          }
        }
        for k in i + 1 .. j {
          let ml_free_energy = ss_max_free_energy_mats.max_free_energy_mat_4_at_least_1_base_pairings_on_mls[i + 1][k - 1] + ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[k][j - 1] + CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop), invert_bp(&(seq[i + 1], seq[j - 1])))] + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
          if ml_free_energy > max_free_energy {max_free_energy = ml_free_energy};
        }
        ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings.insert(pp_closing_loop, max_free_energy);
      }
      max_free_energy = NEG_INFINITY;
      for k in i + 1 .. j + 1 {
        let accessible_pp = (i, k);
        let accessible_bp = (seq[i], seq[k]);
        if !ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
        let ss_max_free_energy_4_bp = ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings[&accessible_pp];
        let free_energy_4_rightmost_base_pairing = ss_max_free_energy_4_bp + (if i > 0 && k < seq_len - 1 {
          ML_TM_DELTA_FES[&(accessible_bp, (seq[i - 1], seq[k + 1]))]
        } else if i > 0 {
          FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])]
        } else if k < seq_len - 1 {
          THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[k + 1])]
        } else {
          0.
        } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
        if free_energy_4_rightmost_base_pairing > max_free_energy {max_free_energy = free_energy_4_rightmost_base_pairing};
      }
      ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[i][j] = max_free_energy;
      max_free_energy = NEG_INFINITY;
      for k in i .. j {
        let ss_max_free_energy_4_rightmost_base_pairings = ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[k][j];
        if ss_max_free_energy_4_rightmost_base_pairings.is_finite() {
          let free_energy = if i == 0 && k == 0 {0.} else {ss_max_free_energy_mats.max_free_energy_mat[i][k - 1]} + ss_max_free_energy_4_rightmost_base_pairings;
        if free_energy > max_free_energy {max_free_energy = free_energy};
        }
      }
      ss_max_free_energy_mats.max_free_energy_mat[i][j] = max_free_energy;
      max_free_energy = ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[i][j] + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
      for k in i + 1 .. j {
        let ss_max_free_energy_4_rightmost_base_pairings = ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[k][j];
        let free_energy_4_at_least_1_base_pairings_on_mls = ss_max_free_energy_4_rightmost_base_pairings + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        if free_energy_4_at_least_1_base_pairings_on_mls > max_free_energy {max_free_energy = free_energy_4_at_least_1_base_pairings_on_mls};
        let free_energy_4_at_least_1_base_pairings_on_mls = ss_max_free_energy_mats.max_free_energy_mat_4_at_least_1_base_pairings_on_mls[i][k - 1] + ss_max_free_energy_4_rightmost_base_pairings + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        if free_energy_4_at_least_1_base_pairings_on_mls > max_free_energy {max_free_energy = free_energy_4_at_least_1_base_pairings_on_mls};
      }
      ss_max_free_energy_mats.max_free_energy_mat_4_at_least_1_base_pairings_on_mls[i][j] = max_free_energy;
    }
  }
  ss_max_free_energy_mats.max_free_energy_mat[0][seq_len - 1]
}

#[inline]
pub fn get_ss_part_func_mats(seq: SeqSlice, seq_len: usize, invert_exp_max_free_energy: FreeEnergy) -> SsPartFuncMats {
  let mut ss_part_func_mats = SsPartFuncMats::new(seq_len);
  for sub_seq_len in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. seq_len + 1 {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let pp_closing_loop = (i, j);
      let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
      let mut sum = 0.;
      if pp_closing_loop.1 - pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL && is_canonical(&bp_closing_loop) {
        sum += get_exp_hl_fe(seq, &pp_closing_loop) * invert_exp_max_free_energy;
        for k in i + 1 .. j - 1 {
          for l in k + 1 .. j {
            if j - l - 1 + k - i - 1 > MAX_2_LOOP_LEN {continue;}
            let accessible_pp = (k, l);
            if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
            let ss_part_func_4_base_pairing = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
            sum += ss_part_func_4_base_pairing * get_exp_2_loop_fe(seq, &pp_closing_loop, &accessible_pp);
          }
        }
        for k in i + 1 .. j {
          sum += ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[i + 1][k - 1] / invert_exp_max_free_energy * ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][j - 1] * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop), invert_bp(&(seq[i + 1], seq[j - 1])))] * if is_au_or_gu(&bp_closing_loop) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
        }
        ss_part_func_mats.part_func_mat_4_base_pairings.insert(pp_closing_loop, sum);
      }
      sum = 0.;
      for k in i + 1 .. j + 1 {
        let accessible_pp = (i, k);
        let accessible_bp = (seq[i], seq[k]);
        if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
        let ss_part_func_4_bp = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
        sum += ss_part_func_4_bp * (if i > 0 && k < seq_len - 1 {
          EXP_ML_TM_DELTA_FES[&(accessible_bp, (seq[i - 1], seq[k + 1]))]
        } else if i > 0 {
          EXP_FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])]
        } else if k < seq_len - 1 {
          EXP_THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[k + 1])]
        } else {
          1.
        } * if is_au_or_gu(&accessible_bp) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.});
      }
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[i][j] = sum;
      sum = 1. * invert_exp_max_free_energy;
      for k in i .. j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][j];
        if ss_part_func_4_rightmost_base_pairings == 0. {
          continue;
        }
        sum += if i == 0 && k == 0 {1.} else {ss_part_func_mats.part_func_flag_mat[i][k - 1].part_func} / if (i == 0 && k == 0) || !ss_part_func_mats.part_func_flag_mat[i][k - 1].flag {1.} else {invert_exp_max_free_energy} * ss_part_func_4_rightmost_base_pairings;
      }
      ss_part_func_mats.part_func_flag_mat[i][j].part_func = sum;
      ss_part_func_mats.part_func_flag_mat[i][j].flag = sum != 0.;
      sum = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[i][j] * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
      for k in i + 1 .. j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][j];
        sum += ss_part_func_4_rightmost_base_pairings * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
        sum += ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[i][k - 1] / invert_exp_max_free_energy * ss_part_func_4_rightmost_base_pairings * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
      }
      ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[i][j] = sum;
    }
  }
  ss_part_func_mats
}

#[inline]
fn get_base_pairing_prob_mat(seq: SeqSlice, ss_part_func_mats: &SsPartFuncMats, seq_len: usize, invert_exp_max_free_energy: FreeEnergy) -> SparseProbMat {
  let ss_part_func = ss_part_func_mats.part_func_flag_mat[0][seq_len - 1].part_func;
  let mut bpp_mat = SparseProbMat::default();
  let mut prob_mat_4_mls_1 = vec![vec![0.; seq_len]; seq_len];
  let mut prob_mat_4_mls_2 = prob_mat_4_mls_1.clone();
  for sub_seq_len in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. seq_len + 1).rev() {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let mut sum_1 = 0.;
      let mut sum_2 = sum_1;
      for k in j + 1 .. seq_len {
        let pp_closing_loop = (i, k);
        if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&pp_closing_loop) {continue;}
        let ss_part_func_4_base_pairing = ss_part_func_mats.part_func_mat_4_base_pairings[&pp_closing_loop];
        let bpp = bpp_mat[&pp_closing_loop];
        let bp_closing_loop = (seq[i], seq[k]);
        let exp_ml_tm_delta_fe = EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop), invert_bp(&(seq[i + 1], seq[k - 1])))];
        let coefficient = bpp * exp_ml_tm_delta_fe * if is_au_or_gu(&bp_closing_loop) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.} / ss_part_func_4_base_pairing;
        sum_1 += coefficient * ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[j + 1][pp_closing_loop.1 - 1] / invert_exp_max_free_energy;
        sum_2 += coefficient;
      }
      prob_mat_4_mls_1[i][j] = sum_1;
      prob_mat_4_mls_2[i][j] = sum_2;
      let accessible_pp = (i, j);
      let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
      if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
      let ss_part_func_4_base_pairing_1 = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
      let mut sum = if accessible_pp.0 < 1 {1.} else {ss_part_func_mats.part_func_flag_mat[0][accessible_pp.0 - 1].part_func} / if accessible_pp.0 < 1 || !ss_part_func_mats.part_func_flag_mat[0][accessible_pp.0 - 1].flag {1.} else {invert_exp_max_free_energy} * ss_part_func_4_base_pairing_1 * if accessible_pp.1 > seq_len - 2 {1.} else {ss_part_func_mats.part_func_flag_mat[accessible_pp.1 + 1][seq_len - 1].part_func} / if accessible_pp.1 > seq_len - 2 || !ss_part_func_mats.part_func_flag_mat[accessible_pp.1 + 1][seq_len - 1].flag {1.} else {invert_exp_max_free_energy} / ss_part_func * (if i > 0 && j < seq_len - 1 {
        EXP_ML_TM_DELTA_FES[&(accessible_bp, (seq[i - 1], seq[j + 1]))]
      } else if i > 0 {
        EXP_FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])]
      } else if j < seq_len - 1 {
        EXP_THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[j + 1])]
      } else {
        1.
      } * if is_au_or_gu(&accessible_bp) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.});
      for k in 0 .. i {
        for l in j + 1 .. seq_len {
          if l - j - 1 + i - k - 1 > MAX_2_LOOP_LEN {continue;}
          let pp_closing_loop = (k, l);
          if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&pp_closing_loop) {continue;}
          let ss_part_func_4_base_pairing_2 = ss_part_func_mats.part_func_mat_4_base_pairings[&pp_closing_loop];
          sum += bpp_mat[&pp_closing_loop] * ss_part_func_4_base_pairing_1 / ss_part_func_4_base_pairing_2 * get_exp_2_loop_fe(seq, &pp_closing_loop, &accessible_pp) as FreeEnergy;
        }
      }
      let coefficient = ss_part_func_4_base_pairing_1 * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * if i > 0 && j < seq_len - 1 {
        EXP_ML_TM_DELTA_FES[&(accessible_bp, (seq[i - 1], seq[j + 1]))]
      } else if i > 0 {
        EXP_FIVE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[i - 1])]
      } else if j < seq_len - 1 {
        EXP_THREE_PRIME_DE_DELTA_FES[&(accessible_bp, seq[j + 1])]
      } else {
        1.
      } * if is_au_or_gu(&accessible_bp) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
      for k in 0 .. i {
        let ss_part_func_4_at_least_1_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][i - 1] / invert_exp_max_free_energy;
        sum += coefficient * prob_mat_4_mls_2[k][j] * ss_part_func_4_at_least_1_base_pairings_on_mls;
        let prob_4_mls = prob_mat_4_mls_1[k][j];
        sum += coefficient * prob_4_mls;
        sum += coefficient * prob_4_mls * ss_part_func_4_at_least_1_base_pairings_on_mls;
      }
      bpp_mat.insert(accessible_pp, sum);
    }
  }
  bpp_mat
}

#[inline]
pub fn logsumexp(sum: &mut FreeEnergy, new_term: FreeEnergy) {
  *sum = if !sum.is_finite() {
   new_term
  } else {
    let max = sum.max(new_term);
    max + ((if *sum == max {new_term - max} else {*sum - max}).exp() + 1.).ln()
  };
}
