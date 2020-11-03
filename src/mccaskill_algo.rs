use utils::*;

pub struct SsPartFuncMats<T: Hash> {
  pub part_func_mat: PartFuncMat,
  pub part_func_mat_4_rightmost_base_pairings: PartFuncMat,
  pub part_func_mat_4_base_pairings: HashMap<PosPair<T>, PartFunc>,
  pub part_func_mat_4_at_least_1_base_pairings_on_mls: PartFuncMat,
}
pub type FreeEnergies = Vec<FreeEnergy>;
pub type FreeEnergyMat = Vec<FreeEnergies>;
pub type SparseFreeEnergyMat<T> = HashMap<PosPair<T>, FreeEnergy>;
pub type SparsePartFuncMat<T> = HashMap<PosPair<T>, PartFunc>;
pub type SparseProbMat<T> = HashMap<PosPair<T>, Prob>;
#[derive(Clone)]
pub struct SsFreeEnergyMats<T: Hash> {
  pub hl_fe_mat: HashMap<PosPair<T>, FreeEnergy>,
  pub twoloop_fe_4d_mat: HashMap<PosQuadruple<T>, FreeEnergy>,
}
pub type FreeEnergy4dMat<T> = HashMap<PosQuadruple<T>, FreeEnergy>;
#[derive(Clone)]
pub struct SsProbMats<T: Hash> {
  pub bpp_mat: SparseProbMat<T>,
  pub access_bpp_mat_4_2l: SparseProbMat<T>,
  pub access_bpp_mat_4_ml: SparseProbMat<T>,
  pub bpp_mat_4_el: SparseProbMat<T>,
}

impl<T: Hash> SsPartFuncMats<T> {
  fn new(seq_len: usize) -> SsPartFuncMats<T> {
    let neg_inf_mat = vec![vec![NEG_INFINITY; seq_len]; seq_len];
    SsPartFuncMats {
      part_func_mat: vec![vec![0.; seq_len]; seq_len],
      part_func_mat_4_rightmost_base_pairings: neg_inf_mat.clone(),
      part_func_mat_4_base_pairings: SparsePartFuncMat::<T>::default(),
      part_func_mat_4_at_least_1_base_pairings_on_mls: neg_inf_mat,
    }
  }
}

impl<T: Unsigned + PrimInt + Hash + One> SsFreeEnergyMats<T> {
  pub fn new() -> SsFreeEnergyMats<T> {
    let free_energy_mat = SparseFreeEnergyMat::<T>::default();
    let free_energy_4d_mat = FreeEnergy4dMat::<T>::default();
    SsFreeEnergyMats {
      hl_fe_mat: free_energy_mat,
      twoloop_fe_4d_mat: free_energy_4d_mat,
    }
  }
  pub fn sparsify(&mut self, bpp_mat: &SparseProbMat<T>, min_bpp: Prob) {
    self.hl_fe_mat = self.hl_fe_mat.iter().filter(|(pos_pair, _)| {bpp_mat[pos_pair] >= min_bpp}).map(|(&(i, j), &free_energy)| {((i + T::one(), j + T::one()), free_energy)}).collect();
    self.twoloop_fe_4d_mat = self.twoloop_fe_4d_mat.iter().filter(|(&(i, j, k, l), _)| {bpp_mat[&(i, j)] >= min_bpp && bpp_mat[&(k, l)] >= min_bpp}).map(|(&(i, j, k, l), &free_energy)| {((i + T::one(), j + T::one(), k + T::one(), l + T::one()), free_energy)}).collect();
  }
}

impl<T: Unsigned + PrimInt + Hash> SsProbMats<T> {
  pub fn new() -> SsProbMats<T> {
    let prob_mat = HashMap::<PosPair<T>, Prob>::default();
    SsProbMats {
      bpp_mat: prob_mat.clone(),
      access_bpp_mat_4_2l: prob_mat.clone(),
      access_bpp_mat_4_ml: prob_mat.clone(),
      bpp_mat_4_el: prob_mat,
    }
  }
}

pub fn mccaskill_algo<T>(seq: SeqSlice) -> (SsProbMats<T>, SsFreeEnergyMats<T>)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let seq_len = seq.len();
  let mut ss_free_energy_mats = SsFreeEnergyMats::<T>::new();
  let ss_part_func_mats = get_ss_part_func_mats::<T>(seq, seq_len, &mut ss_free_energy_mats);
  let bpp_mats = get_base_pairing_prob_mats::<T>(seq, &ss_part_func_mats, seq_len, &ss_free_energy_mats);
  (bpp_mats, ss_free_energy_mats)
}

pub fn get_bpp_and_unpair_prob_mats<T>(seq: SeqSlice) -> (SsProbMats<T>, Probs, SsFreeEnergyMats<T>)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let seq_len = seq.len();
  let (bpp_mats, ss_free_energy_mats) = mccaskill_algo::<T>(&seq[..]);
  let mut unpair_prob_mat = vec![1.; seq_len];
  let seq_len = T::from_usize(seq_len).unwrap();
  for i in range(T::zero(), seq_len) {
    let long_i = i.to_usize().unwrap();
    let mut sum = 0.;
    for j in range(T::zero(), seq_len) {
      if j == i {continue;}
      let pp = if j < i {(j, i)} else {(i, j)};
      if !bpp_mats.bpp_mat.contains_key(&pp) {continue;}
      sum += bpp_mats.bpp_mat[&pp];
    }
    debug_assert!(0. <= sum && sum <= 1.);
    unpair_prob_mat[long_i] = 1. - sum;
  }
  (bpp_mats, unpair_prob_mat, ss_free_energy_mats)
}

pub fn get_ss_part_func_mats<T>(seq: SeqSlice, seq_len: usize, ss_free_energy_mats: &mut SsFreeEnergyMats<T>) -> SsPartFuncMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let mut ss_part_func_mats = SsPartFuncMats::<T>::new(seq_len);
  let seq_len = T::from_usize(seq_len).unwrap();
  for sub_seq_len in range_inclusive(T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(), seq_len) {
    for i in range_inclusive(T::zero(), seq_len - sub_seq_len) {
      let j = i + sub_seq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let pp_closing_loop = (i, j);
      let long_pp_closing_loop = (long_i, long_j);
      let bp_closing_loop = (seq[long_i], seq[long_j]);
      let mut sum = NEG_INFINITY;
      if long_pp_closing_loop.1 - long_pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL && is_canonical(&bp_closing_loop) {
        let hl_fe = get_hl_fe(seq, &long_pp_closing_loop);
        ss_free_energy_mats.hl_fe_mat.insert(pp_closing_loop, hl_fe);
        logsumexp(&mut sum, hl_fe);
        for k in range(i + T::one(), j - T::one()) {
          let long_k = k.to_usize().unwrap();
          for l in range(k + T::one(), j) {
            let long_l = l.to_usize().unwrap();
            if long_j - long_l - 1 + long_k - long_i - 1 > MAX_2_LOOP_LEN {continue;}
            let accessible_pp = (k, l);
            let long_accessible_pp = (long_k, long_l);
            if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
            let ss_part_func_4_base_pairing = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
            let twoloop_free_energy = get_2_loop_fe(seq, &long_pp_closing_loop, &long_accessible_pp);
            ss_free_energy_mats.twoloop_fe_4d_mat.insert((i, j, k, l), twoloop_free_energy);
            logsumexp(&mut sum, ss_part_func_4_base_pairing + twoloop_free_energy);
          }
        }
        let invert_bp_closing_loop = invert_bp(&bp_closing_loop);
        let invert_stacking_bp = invert_bp(&(seq[long_i + 1], seq[long_j - 1]));
        let ml_tm_delta_fe = ML_TM_DELTA_FES[invert_bp_closing_loop.0][invert_bp_closing_loop.1][invert_stacking_bp.0][invert_stacking_bp.1];
        for k in long_i + 1 .. long_j {
          logsumexp(&mut sum, ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i + 1][k - 1] + ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j - 1] + CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ml_tm_delta_fe + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
        }
        ss_part_func_mats.part_func_mat_4_base_pairings.insert(pp_closing_loop, sum);
      }
      sum = NEG_INFINITY;
      for k in range_inclusive(i + T::one(), j) {
        let long_k = k.to_usize().unwrap();
        let accessible_pp = (i, k);
        let accessible_bp = (seq[long_i], seq[long_k]);
        if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
        let ss_part_func_4_bp = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
        logsumexp(&mut sum, ss_part_func_4_bp + if i > T::zero() && k < seq_len - T::one() {
          ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]][seq[long_k + 1]]
        } else if i > T::zero() {
          FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]]
        } else if k < seq_len - T::one() {
          THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_k + 1]]
        } else {
          0.
        } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.});
      }
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[long_i][long_j] = sum;
      sum = 0.;
      for k in long_i .. long_j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j];
        if ss_part_func_4_rightmost_base_pairings == NEG_INFINITY {
          continue;
        }
        let part_func = if long_i == 0 && k == 0 {0.} else {ss_part_func_mats.part_func_mat[long_i][k - 1]};
        logsumexp(&mut sum, part_func + ss_part_func_4_rightmost_base_pairings);
      }
      ss_part_func_mats.part_func_mat[long_i][long_j] = sum;
      sum = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[long_i][long_j] + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
      for k in long_i + 1 .. long_j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j];
        logsumexp(&mut sum, ss_part_func_4_rightmost_base_pairings + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
        logsumexp(&mut sum, ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][k - 1] + ss_part_func_4_rightmost_base_pairings + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
      }
      ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][long_j] = sum;
    }
  }
  ss_part_func_mats
}

fn get_base_pairing_prob_mats<T>(seq: SeqSlice, ss_part_func_mats: &SsPartFuncMats<T>, seq_len: usize, ss_free_energy_mats: &SsFreeEnergyMats<T>) -> SsProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let ss_part_func = ss_part_func_mats.part_func_mat[0][seq_len - 1];
  let mut bpp_mats = SsProbMats::<T>::new();
  let mut prob_mat_4_mls_1 = vec![vec![NEG_INFINITY; seq_len]; seq_len];
  let mut prob_mat_4_mls_2 = prob_mat_4_mls_1.clone();
  let short_seq_len = T::from_usize(seq_len).unwrap();
  for sub_seq_len in range_inclusive(T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(), short_seq_len).rev() {
    for i in range_inclusive(T::zero(), short_seq_len - sub_seq_len) {
      let j = i + sub_seq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let mut sum_1 = NEG_INFINITY;
      let mut sum_2 = sum_1;
      for k in range(j + T::one(), short_seq_len) {
        let long_k = k.to_usize().unwrap();
        let pp_closing_loop = (i, k);
        if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&pp_closing_loop) {continue;}
        let ss_part_func_4_base_pairing = ss_part_func_mats.part_func_mat_4_base_pairings[&pp_closing_loop];
        let bpp = bpp_mats.bpp_mat[&pp_closing_loop];
        let bp_closing_loop = (seq[long_i], seq[long_k]);
        let invert_bp_closing_loop = invert_bp(&bp_closing_loop);
        let invert_stacking_bp = invert_bp(&(seq[long_i + 1], seq[long_k - 1]));
        let ml_tm_delta_fe = ML_TM_DELTA_FES[invert_bp_closing_loop.0][invert_bp_closing_loop.1][invert_stacking_bp.0][invert_stacking_bp.1];
        let coefficient = bpp + ml_tm_delta_fe + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.} - ss_part_func_4_base_pairing;
        logsumexp(&mut sum_1, coefficient + ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_j + 1][long_k - 1]);
        logsumexp(&mut sum_2, coefficient);
      }
      prob_mat_4_mls_1[long_i][long_j] = sum_1;
      prob_mat_4_mls_2[long_i][long_j] = sum_2;
      let accessible_pp = (i, j);
      let accessible_bp = (seq[long_i], seq[long_j]);
      if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&accessible_pp) {continue;}
      let ss_part_func_4_base_pairing_1 = ss_part_func_mats.part_func_mat_4_base_pairings[&accessible_pp];
      let part_func_pair = (
        if accessible_pp.0 < T::one() {0.} else {ss_part_func_mats.part_func_mat[0][long_i - 1]},
        if accessible_pp.1 > short_seq_len - T::from_usize(2).unwrap() {0.} else {ss_part_func_mats.part_func_mat[long_j + 1][seq_len - 1]},
      );
      let mut sum = part_func_pair.0 + ss_part_func_4_base_pairing_1 + part_func_pair.1 - ss_part_func + if i > T::zero() && j < short_seq_len - T::one() {
        ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]][seq[long_j + 1]]
      } else if i > T::zero() {
        FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]]
      } else if j < short_seq_len - T::one() {
        THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_j + 1]]
      } else {
        0.
      } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
      if sum > NEG_INFINITY {
        bpp_mats.bpp_mat_4_el.insert(accessible_pp, sum.exp());
      }
      let mut bpp_4_2l = NEG_INFINITY;
      for k in range(T::zero(), i) {
        let long_k = k.to_usize().unwrap();
        for l in range(j + T::one(), short_seq_len) {
          let long_l = l.to_usize().unwrap();
          if long_l - long_j - 1 + long_i - long_k - 1 > MAX_2_LOOP_LEN {continue;}
          let pp_closing_loop = (k, l);
          if !ss_part_func_mats.part_func_mat_4_base_pairings.contains_key(&pp_closing_loop) {continue;}
          let ss_part_func_4_base_pairing_2 = ss_part_func_mats.part_func_mat_4_base_pairings[&pp_closing_loop];
          logsumexp(&mut bpp_4_2l, bpp_mats.bpp_mat[&pp_closing_loop] + ss_part_func_4_base_pairing_1 - ss_part_func_4_base_pairing_2 + ss_free_energy_mats.twoloop_fe_4d_mat[&(k, l, i, j)]);
        }
      }
      if bpp_4_2l > NEG_INFINITY {
        bpp_mats.access_bpp_mat_4_2l.insert(accessible_pp, bpp_4_2l.exp());
        logsumexp(&mut sum, bpp_4_2l);
      }
      let coefficient = ss_part_func_4_base_pairing_1 + CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + if i > T::zero() && j < short_seq_len - T::one() {
        ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]][seq[long_j + 1]]
      } else if i > T::zero() {
        FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_i - 1]]
      } else if j < short_seq_len - T::one() {
        THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[long_j + 1]]
      } else {
        0.
      } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
      let mut bpp_4_ml = NEG_INFINITY;
      for k in 0 .. long_i {
        let ss_part_func_4_at_least_1_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_i - 1];
        logsumexp(&mut bpp_4_ml, coefficient + prob_mat_4_mls_2[k][long_j] + ss_part_func_4_at_least_1_base_pairings_on_mls);
        let prob_4_mls = prob_mat_4_mls_1[k][long_j];
        logsumexp(&mut bpp_4_ml, coefficient + prob_4_mls);
        logsumexp(&mut bpp_4_ml, coefficient + prob_4_mls + ss_part_func_4_at_least_1_base_pairings_on_mls);
      }
      if bpp_4_ml > NEG_INFINITY {
        bpp_mats.access_bpp_mat_4_ml.insert(accessible_pp, bpp_4_ml.exp());
        logsumexp(&mut sum, bpp_4_ml);
      }
      debug_assert!(NEG_INFINITY <= sum && sum <= 0.);
      bpp_mats.bpp_mat.insert(accessible_pp, sum);
    }
  }
  bpp_mats.bpp_mat = bpp_mats.bpp_mat.iter().map(|(pos_pair, &bpp)| {(*pos_pair, bpp.exp())}).collect();
  bpp_mats
}

pub fn logsumexp(sum: &mut FreeEnergy, new_term: FreeEnergy) {
  *sum = if !sum.is_finite() {
   new_term
  } else {
    let max = sum.max(new_term);
    max + ((if *sum == max {new_term - max} else {*sum - max}).exp() + 1.).ln()
  };
}
