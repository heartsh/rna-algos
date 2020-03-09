use utils::*;

pub struct SsPartFuncMats {
  pub part_func_mat: PartFuncMat,
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
pub type SparseFreeEnergyMat = FxHashMap<PosPair, FreeEnergy>;
pub type SparsePartFuncMat = FxHashMap<PosPair, PartFunc>;
pub type SparseProbMat = FxHashMap<PosPair, Prob>;
#[derive(Clone)]
pub struct SsFreeEnergyMats {
  pub hl_fe_mat: SparseFreeEnergyMat,
  pub exp_hl_fe_mat: SparseFreeEnergyMat,
  pub twoloop_fe_4d_mat: FreeEnergy4dMat,
  pub exp_2loop_fe_4d_mat: FreeEnergy4dMat,
}
pub type PosQuadruple = (Pos, Pos, Pos, Pos);
pub type FreeEnergy4dMat = FxHashMap<PosQuadruple, FreeEnergy>;

impl SsPartFuncMats {
  fn new(seq_len: usize) -> SsPartFuncMats {
    let zero_mat = vec![vec![0.; seq_len]; seq_len];
    SsPartFuncMats {
      part_func_mat: vec![vec![1.; seq_len]; seq_len],
      part_func_mat_4_rightmost_base_pairings: zero_mat.clone(),
      part_func_mat_4_base_pairings: SparsePartFuncMat::default(),
      part_func_mat_4_at_least_1_base_pairings_on_mls: zero_mat,
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

impl SsFreeEnergyMats {
  pub fn new() -> SsFreeEnergyMats {
    let free_energy_mat = SparseFreeEnergyMat::default();
    let free_energy_4d_mat = FreeEnergy4dMat::default();
    SsFreeEnergyMats {
      hl_fe_mat: free_energy_mat.clone(),
      exp_hl_fe_mat: free_energy_mat,
      twoloop_fe_4d_mat: free_energy_4d_mat.clone(),
      exp_2loop_fe_4d_mat: free_energy_4d_mat,
    }
  }
  pub fn sparsify(&mut self, bpp_mat: &SparseProbMat, min_bpp: Prob) {
    self.hl_fe_mat = self.hl_fe_mat.iter().filter(|(pos_pair, _)| {bpp_mat[pos_pair] >= min_bpp}).map(|(&(i, j), &free_energy)| {((i + 1, j + 1), free_energy)}).collect();
    self.exp_hl_fe_mat = self.exp_hl_fe_mat.iter().filter(|(pos_pair, _)| {bpp_mat[pos_pair] >= min_bpp}).map(|(&(i, j), &free_energy)| {((i + 1, j + 1), free_energy)}).collect();
    self.twoloop_fe_4d_mat = self.twoloop_fe_4d_mat.iter().filter(|(&(i, j, k, l), _)| {bpp_mat[&(i, j)] >= min_bpp && bpp_mat[&(k, l)] >= min_bpp}).map(|(&(i, j, k, l), &free_energy)| {((i + 1, j + 1, k + 1, l + 1), free_energy)}).collect();
    self.exp_2loop_fe_4d_mat = self.exp_2loop_fe_4d_mat.iter().filter(|(&(i, j, k, l), _)| {bpp_mat[&(i, j)] >= min_bpp && bpp_mat[&(k, l)] >= min_bpp}).map(|(&(i, j, k, l), &free_energy)| {((i + 1, j + 1, k + 1, l + 1), free_energy)}).collect();
  }
}

pub const CONST_4_INIT_ML_DELTA_FE: FreeEnergy = - INVERSE_TEMPERATURE * 9.3;
pub const COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE: FreeEnergy = - INVERSE_TEMPERATURE * (-0.9);
lazy_static! {
  pub static ref EXP_CONST_4_INIT_ML_DELTA_FE: FreeEnergy = CONST_4_INIT_ML_DELTA_FE.exp();
  pub static ref EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE: FreeEnergy = COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE.exp();
}

pub fn mccaskill_algo(seq: SeqSlice) -> (SparseProbMat, FreeEnergy, SsFreeEnergyMats) {
  let seq_len = seq.len();
  let mut ss_free_energy_mats = SsFreeEnergyMats::new();
  let max_free_energy = get_max_free_energy(seq, seq_len, &mut ss_free_energy_mats);
  let invert_exp_max_free_energy = 1. / max_free_energy.exp();
  let ss_part_func_mats = get_ss_part_func_mats(seq, seq_len, invert_exp_max_free_energy, &mut ss_free_energy_mats);
  let bpp_mat = get_base_pairing_prob_mat(seq, &ss_part_func_mats, seq_len, invert_exp_max_free_energy);
  (bpp_mat, max_free_energy, ss_free_energy_mats)
}

pub fn get_bpp_and_unpair_prob_mats(seq: SeqSlice) -> (SparseProbMat, Probs, FreeEnergy, SsFreeEnergyMats) {
  let seq_len = seq.len();
  let (bpp_mat, max_free_energy, ss_free_energy_mats) = mccaskill_algo(&seq[..]);
  let mut unpair_prob_mat = vec![1.; seq_len];
  let seq_len = seq_len as Pos;
  for i in 0 .. seq_len {
    let long_i = i as usize;
    let mut sum = 0.;
    for j in 0 .. seq_len {
      if j == i {continue;}
      let pp = if j < i {(j, i)} else {(i, j)};
      if !bpp_mat.contains_key(&pp) {continue;}
      sum += bpp_mat[&pp];
    }
    assert!(0. <= sum && sum <= 1.);
    unpair_prob_mat[long_i] = 1. - sum;
  }
  (bpp_mat, unpair_prob_mat, max_free_energy, ss_free_energy_mats)
}

pub fn get_max_free_energy(seq: SeqSlice, seq_len: usize, ss_free_energy_mats: &mut SsFreeEnergyMats) -> FreeEnergy {
  let mut ss_max_free_energy_mats = SsMaxFreeEnergyMats::new(seq_len);
  let seq_len = seq_len as Pos;
  for sub_seq_len in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. seq_len + 1 {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let (long_i, long_j) = (i as usize, j as usize);
      let pp_closing_loop = (i, j);
      let long_pp_closing_loop = (long_i, long_j);
      let bp_closing_loop = (seq[long_i], seq[long_j]);
      let mut max_free_energy;
      if pp_closing_loop.1 - pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos && is_canonical(&bp_closing_loop) {
        max_free_energy = get_hl_fe(seq, &long_pp_closing_loop);
        ss_free_energy_mats.hl_fe_mat.insert(pp_closing_loop, max_free_energy);
        for k in i + 1 .. j - 1 {
          let long_k = k as usize;
          for l in k + 1 .. j {
            let long_l = l as usize;
            if long_j - long_l - 1 + long_k - long_i - 1 > MAX_2_LOOP_LEN {continue;}
            let accessible_pp = (k, l);
            let long_accessible_pp = (long_k, long_l);
            if !ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
            let ss_max_free_energy_4_base_pairing = ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings[&accessible_pp];
            // let twoloop_free_energy = ss_max_free_energy_4_base_pairing + get_2_loop_fe(seq, &long_pp_closing_loop, &long_accessible_pp);
            let twoloop_free_energy = get_2_loop_fe(seq, &long_pp_closing_loop, &long_accessible_pp);
            ss_free_energy_mats.twoloop_fe_4d_mat.insert((i, j, k, l), twoloop_free_energy);
            let twoloop_free_energy = ss_max_free_energy_4_base_pairing + twoloop_free_energy;
            if twoloop_free_energy > max_free_energy {max_free_energy = twoloop_free_energy};
          }
        }
        let invert_bp_closing_loop = invert_bp(&bp_closing_loop);
        let invert_stacking_bp = invert_bp(&(seq[long_i + 1], seq[long_j - 1]));
        for k in long_i + 1 .. long_j {
          let ml_free_energy = ss_max_free_energy_mats.max_free_energy_mat_4_at_least_1_base_pairings_on_mls[long_i + 1][k - 1] + ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[k][long_j - 1] + CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ML_TM_DELTA_FES[invert_bp_closing_loop.0][invert_bp_closing_loop.1][invert_stacking_bp.0][invert_stacking_bp.1] + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
          if ml_free_energy > max_free_energy {max_free_energy = ml_free_energy};
        }
        ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings.insert(pp_closing_loop, max_free_energy);
      }
      max_free_energy = NEG_INFINITY;
      for k in i + 1 .. j + 1 {
        let long_k = k as usize;
        let accessible_pp = (i, k);
        let accessible_bp = (seq[long_i], seq[long_k]);
        if !ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
        let ss_max_free_energy_4_bp = ss_max_free_energy_mats.max_free_energy_mat_4_base_pairings[&accessible_pp];
        let free_energy_4_rightmost_base_pairing = ss_max_free_energy_4_bp + (if i > 0 && k < seq_len - 1 {
          ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]][seq[long_k + 1]]
        } else if i > 0 {
          FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]]
        } else if k < seq_len - 1 {
          THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_k + 1]]
        } else {
          0.
        } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
        if free_energy_4_rightmost_base_pairing > max_free_energy {max_free_energy = free_energy_4_rightmost_base_pairing};
      }
      ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[long_i][long_j] = max_free_energy;
      max_free_energy = 0.;
      for k in long_i .. long_j {
        let ss_max_free_energy_4_rightmost_base_pairings = ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[k][long_j];
        if ss_max_free_energy_4_rightmost_base_pairings.is_finite() {
          let free_energy = if i == 0 && k == 0 {0.} else {ss_max_free_energy_mats.max_free_energy_mat[long_i][k - 1]} + ss_max_free_energy_4_rightmost_base_pairings;
        if free_energy > max_free_energy {max_free_energy = free_energy};
        }
      }
      ss_max_free_energy_mats.max_free_energy_mat[long_i][long_j] = max_free_energy;
      max_free_energy = ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[long_i][long_j] + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
      for k in long_i + 1 .. long_j {
        let ss_max_free_energy_4_rightmost_base_pairings = ss_max_free_energy_mats.max_free_energy_mat_4_rightmost_base_pairings[k][long_j];
        let free_energy_4_at_least_1_base_pairings_on_mls = ss_max_free_energy_4_rightmost_base_pairings + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        if free_energy_4_at_least_1_base_pairings_on_mls > max_free_energy {max_free_energy = free_energy_4_at_least_1_base_pairings_on_mls};
        let free_energy_4_at_least_1_base_pairings_on_mls = ss_max_free_energy_mats.max_free_energy_mat_4_at_least_1_base_pairings_on_mls[long_i][k - 1] + ss_max_free_energy_4_rightmost_base_pairings + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        if free_energy_4_at_least_1_base_pairings_on_mls > max_free_energy {max_free_energy = free_energy_4_at_least_1_base_pairings_on_mls};
      }
      ss_max_free_energy_mats.max_free_energy_mat_4_at_least_1_base_pairings_on_mls[long_i][long_j] = max_free_energy;
    }
  }
  ss_max_free_energy_mats.max_free_energy_mat[0][seq_len as usize - 1]
}

pub fn get_ss_part_func_mats(seq: SeqSlice, seq_len: usize, invert_exp_max_free_energy: FreeEnergy, ss_free_energy_mats: &mut SsFreeEnergyMats) -> SsPartFuncMats {
  let mut ss_part_func_mats = SsPartFuncMats::new(seq_len);
  let seq_len = seq_len as Pos;
  for sub_seq_len in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. seq_len + 1 {
    for i in 0 .. seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let (long_i, long_j) = (i as usize, j as usize);
      let pp_closing_loop = (i, j);
      let long_pp_closing_loop = (long_i, long_j);
      let bp_closing_loop = (seq[long_i], seq[long_j]);
      let mut sum = 0.;
      if long_pp_closing_loop.1 - long_pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL && is_canonical(&bp_closing_loop) {
        let exp_hl_fe = get_exp_hl_fe(seq, &long_pp_closing_loop);
        ss_free_energy_mats.exp_hl_fe_mat.insert(pp_closing_loop, exp_hl_fe);
        sum += exp_hl_fe * invert_exp_max_free_energy;
        for k in i + 1 .. j - 1 {
          let long_k = k as usize;
          for l in k + 1 .. j {
            let long_l = l as usize;
            if long_j - long_l - 1 + long_k - long_i - 1 > MAX_2_LOOP_LEN {continue;}
            let accessible_pp = (k, l);
            let long_accessible_pp = (long_k, long_l);
            if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
            let ss_part_func_4_base_pairing = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
            let exp_2loop_free_energy = get_exp_2_loop_fe(seq, &long_pp_closing_loop, &long_accessible_pp);
            ss_free_energy_mats.exp_2loop_fe_4d_mat.insert((i, j, k, l), exp_2loop_free_energy);
            sum += ss_part_func_4_base_pairing * exp_2loop_free_energy;
          }
        }
        let invert_bp_closing_loop = invert_bp(&bp_closing_loop);
        let invert_stacking_bp = invert_bp(&(seq[long_i + 1], seq[long_j - 1]));
        let exp_ml_tm_delta_fe = EXP_ML_TM_DELTA_FES[invert_bp_closing_loop.0][invert_bp_closing_loop.1][invert_stacking_bp.0][invert_stacking_bp.1];
        for k in long_i + 1 .. long_j {
          sum += ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i + 1][k - 1] * ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j - 1] / invert_exp_max_free_energy * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe * if is_au_or_gu(&bp_closing_loop) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
        }
        ss_part_func_mats.part_func_mat_4_base_pairings.insert(pp_closing_loop, sum);
      }
      sum = 0.;
      for k in i + 1 .. j + 1 {
        let long_k = k as usize;
        let accessible_pp = (i, k);
        let accessible_bp = (seq[long_i], seq[long_k]);
        if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
        let ss_part_func_4_bp = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
        sum += ss_part_func_4_bp * (if i > 0 && k < seq_len - 1 {
          EXP_ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]][seq[long_k + 1]]
        } else if i > 0 {
          EXP_FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]]
        } else if k < seq_len - 1 {
          EXP_THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_k + 1]]
        } else {
          1.
        } * if is_au_or_gu(&accessible_bp) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.});
      }
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[long_i][long_j] = sum;
      sum = invert_exp_max_free_energy;
      for k in long_i .. long_j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j];
        if ss_part_func_4_rightmost_base_pairings == 0. {
          continue;
        }
        let part_func = if i == 0 && k == 0 {1.} else {ss_part_func_mats.part_func_mat[long_i][k - 1]};
        sum += part_func * ss_part_func_4_rightmost_base_pairings / if part_func == 1. {1.} else {invert_exp_max_free_energy};
      }
      ss_part_func_mats.part_func_mat[long_i][long_j] = sum;
      sum = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[long_i][long_j] * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
      for k in long_i + 1 .. long_j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j];
        sum += ss_part_func_4_rightmost_base_pairings * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
        sum += ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][k - 1] * ss_part_func_4_rightmost_base_pairings / invert_exp_max_free_energy * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
      }
      ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][long_j] = sum;
    }
  }
  ss_part_func_mats
}

fn get_base_pairing_prob_mat(seq: SeqSlice, ss_part_func_mats: &SsPartFuncMats, seq_len: usize, invert_exp_max_free_energy: FreeEnergy) -> SparseProbMat {
  let ss_part_func = ss_part_func_mats.part_func_mat[0][seq_len - 1];
  let mut bpp_mat = SparseProbMat::default();
  let mut prob_mat_4_mls_1 = vec![vec![0.; seq_len]; seq_len];
  let mut prob_mat_4_mls_2 = prob_mat_4_mls_1.clone();
  let short_seq_len = seq_len as Pos;
  for sub_seq_len in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. short_seq_len + 1).rev() {
    for i in 0 .. short_seq_len - sub_seq_len + 1 {
      let j = i + sub_seq_len - 1;
      let (long_i, long_j) = (i as usize, j as usize);
      let mut sum_1 = 0.;
      let mut sum_2 = sum_1;
      for k in j + 1 .. short_seq_len {
        let long_k = k as usize;
        let pp_closing_loop = (i, k);
        if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&pp_closing_loop) {continue;}
        let ss_part_func_4_base_pairing = ss_part_func_mats.part_func_mat_4_base_pairings[&pp_closing_loop];
        let bpp = bpp_mat[&pp_closing_loop];
        let bp_closing_loop = (seq[long_i], seq[long_k]);
        let invert_bp_closing_loop = invert_bp(&bp_closing_loop);
        let invert_stacking_bp = invert_bp(&(seq[long_i + 1], seq[long_k - 1]));
        let exp_ml_tm_delta_fe = EXP_ML_TM_DELTA_FES[invert_bp_closing_loop.0][invert_bp_closing_loop.1][invert_stacking_bp.0][invert_stacking_bp.1];
        let coefficient = bpp * exp_ml_tm_delta_fe * if is_au_or_gu(&bp_closing_loop) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.} / ss_part_func_4_base_pairing;
        sum_1 += coefficient * ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_j + 1][long_k - 1];
        sum_2 += coefficient;
      }
      prob_mat_4_mls_1[long_i][long_j] = sum_1;
      prob_mat_4_mls_2[long_i][long_j] = sum_2;
      let accessible_pp = (i, j);
      let long_accessible_pp = (long_i, long_j);
      let accessible_bp = (seq[long_i], seq[long_j]);
      if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
      let ss_part_func_4_base_pairing_1 = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
      let part_func_pair = (
        if accessible_pp.0 < 1 {1.} else {ss_part_func_mats.part_func_mat[0][long_i - 1]},
        if accessible_pp.1 > short_seq_len - 2 {1.} else {ss_part_func_mats.part_func_mat[long_j + 1][seq_len - 1]},
      );
      let mut sum = part_func_pair.0 * ss_part_func_4_base_pairing_1 / if part_func_pair.0 == 1. {1.} else {invert_exp_max_free_energy} * part_func_pair.1 / if part_func_pair.1 == 1. {1.} else {invert_exp_max_free_energy} / ss_part_func * (if i > 0 && j < short_seq_len - 1 {
        EXP_ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]][seq[long_j + 1]]
      } else if i > 0 {
        EXP_FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]]
      } else if j < short_seq_len - 1 {
        EXP_THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_j + 1]]
      } else {
        1.
      } * if is_au_or_gu(&accessible_bp) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.});
      for k in 0 .. i {
        let long_k = k as usize;
        for l in j + 1 .. short_seq_len {
          let long_l = l as usize;
          if long_l - long_j - 1 + long_i - long_k - 1 > MAX_2_LOOP_LEN {continue;}
          let pp_closing_loop = (k, l);
          let long_pp_closing_loop = (long_k, long_l);
          if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&pp_closing_loop) {continue;}
          let ss_part_func_4_base_pairing_2 = ss_part_func_mats.part_func_mat_4_base_pairings[&pp_closing_loop];
          sum += bpp_mat[&pp_closing_loop] * ss_part_func_4_base_pairing_1 / ss_part_func_4_base_pairing_2 * get_exp_2_loop_fe(seq, &long_pp_closing_loop, &long_accessible_pp) as FreeEnergy;
        }
      }
      let coefficient = ss_part_func_4_base_pairing_1 * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * if i > 0 && j < short_seq_len - 1 {
        EXP_ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]][seq[long_j + 1]]
      } else if i > 0 {
        EXP_FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]]
      } else if j < short_seq_len - 1 {
        EXP_THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_j + 1]]
      } else {
        1.
      } * if is_au_or_gu(&accessible_bp) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
      for k in 0 .. long_i {
        let ss_part_func_4_at_least_1_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_i - 1];
        sum += coefficient * prob_mat_4_mls_2[k][long_j] * ss_part_func_4_at_least_1_base_pairings_on_mls / invert_exp_max_free_energy;
        let prob_4_mls = prob_mat_4_mls_1[k][long_j];
        sum += coefficient * prob_4_mls / invert_exp_max_free_energy;
        sum += coefficient * prob_4_mls / invert_exp_max_free_energy * ss_part_func_4_at_least_1_base_pairings_on_mls / invert_exp_max_free_energy;
      }
      assert!(0. <= sum && sum <= 1.);
      bpp_mat.insert(accessible_pp, sum);
    }
  }
  bpp_mat
}

pub fn logsumexp(sum: &mut FreeEnergy, new_term: FreeEnergy) {
  *sum = if !sum.is_finite() {
   new_term
  } else {
    let max = sum.max(new_term);
    max + ((if *sum == max {new_term - max} else {*sum - max}).exp() + 1.).ln()
  };
}
