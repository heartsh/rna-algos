use utils::*;

pub struct SsPartFuncMats<T: Hash> {
  pub part_func_mat: PartFuncMat,
  pub part_func_mat_4_rightmost_base_pairings: PartFuncMat,
  pub part_func_mat_4_rightmost_base_pairings_on_mls: PartFuncMat,
  pub part_func_mat_4_base_pairings: SparsePartFuncMat<T>,
  pub part_func_mat_4_base_pairings_accessible: SparsePartFuncMat<T>,
  pub part_func_mat_4_ml: PartFuncMat,
  pub part_func_mat_4_at_least_1_base_pairings_on_mls: PartFuncMat,
}

#[derive(Clone)]
pub struct SsFreeEnergyMats<T: Hash> {
  pub hl_fe_mat: HashMap<PosPair<T>, FreeEnergy>,
  pub twoloop_fe_4d_mat: HashMap<PosQuadruple<T>, FreeEnergy>,
  pub ml_closing_bp_fe_mat: HashMap<PosPair<T>, FreeEnergy>,
  pub accessible_bp_fe_mat: HashMap<PosPair<T>, FreeEnergy>,
}
pub type FreeEnergy4dMat<T> = HashMap<PosQuadruple<T>, FreeEnergy>;

impl StructFeatureCountSets {
  pub fn new(init_val: FeatureCount) -> StructFeatureCountSets {
    let init_vals = [init_val; NUM_OF_BASES];
    let twod_mat = [init_vals; NUM_OF_BASES];
    let threed_mat = [twod_mat; NUM_OF_BASES];
    let fourd_mat = [threed_mat; NUM_OF_BASES];
    StructFeatureCountSets {
      // The CONTRAfold model.
      hairpin_loop_length_counts: [init_val; CONSPROB_MAX_HAIRPIN_LOOP_LEN + 1],
      bulge_loop_length_counts: [init_val; CONSPROB_MAX_TWOLOOP_LEN],
      interior_loop_length_counts: [init_val; CONSPROB_MAX_TWOLOOP_LEN - 1],
      interior_loop_length_counts_symm: [init_val; CONSPROB_MAX_INTERIOR_LOOP_LEN_SYMM],
      interior_loop_length_counts_asymm: [init_val; CONSPROB_MAX_INTERIOR_LOOP_LEN_ASYMM],
      stack_count_mat: fourd_mat,
      terminal_mismatch_count_mat: fourd_mat,
      left_dangle_count_mat: threed_mat,
      right_dangle_count_mat: threed_mat,
      helix_end_count_mat: twod_mat,
      base_pair_count_mat: twod_mat,
      interior_loop_length_count_mat_explicit: [[init_val; CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT]; CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT],
      bulge_loop_0x1_length_counts: [init_val; NUM_OF_BASES],
      interior_loop_1x1_length_count_mat: [[init_val; NUM_OF_BASES]; NUM_OF_BASES],
      multi_loop_base_count: init_val,
      multi_loop_basepairing_count: init_val,
      multi_loop_accessible_baseunpairing_count: init_val,
      external_loop_accessible_basepairing_count: init_val,
      external_loop_accessible_baseunpairing_count: init_val,
      // The cumulative parameters of the CONTRAfold model.
      hairpin_loop_length_counts_cumulative: [init_val; CONSPROB_MAX_HAIRPIN_LOOP_LEN + 1],
      bulge_loop_length_counts_cumulative: [init_val; CONSPROB_MAX_TWOLOOP_LEN],
      interior_loop_length_counts_cumulative: [init_val; CONSPROB_MAX_TWOLOOP_LEN - 1],
      interior_loop_length_counts_symm_cumulative: [init_val; CONSPROB_MAX_INTERIOR_LOOP_LEN_SYMM],
      interior_loop_length_counts_asymm_cumulative: [init_val; CONSPROB_MAX_INTERIOR_LOOP_LEN_ASYMM],
    }
  }

  pub fn accumulate(&mut self)
  {
    let mut sum = 0.;
    for i in 0 .. self.hairpin_loop_length_counts_cumulative.len() {
      sum += self.hairpin_loop_length_counts[i];
      self.hairpin_loop_length_counts_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0 .. self.bulge_loop_length_counts_cumulative.len() {
      sum += self.bulge_loop_length_counts[i];
      self.bulge_loop_length_counts_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0 .. self.interior_loop_length_counts_cumulative.len() {
      sum += self.interior_loop_length_counts[i];
      self.interior_loop_length_counts_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0 .. self.interior_loop_length_counts_symm_cumulative.len() {
      sum += self.interior_loop_length_counts_symm[i];
      self.interior_loop_length_counts_symm_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0 .. self.interior_loop_length_counts_asymm_cumulative.len() {
      sum += self.interior_loop_length_counts_asymm[i];
      self.interior_loop_length_counts_asymm_cumulative[i] = sum;
    }
  }

  pub fn transfer(&mut self) {
    for (v, &w) in self.hairpin_loop_length_counts.iter_mut().zip(CONTRA_HL_LENGTH_FES_AT_LEAST.iter()) {
      *v = w;
    }
    for (v, &w) in self.bulge_loop_length_counts.iter_mut().zip(CONTRA_BL_LENGTH_FES_AT_LEAST.iter()) {
      *v = w;
    }
    for (v, &w) in self.interior_loop_length_counts.iter_mut().zip(CONTRA_IL_LENGTH_FES_AT_LEAST.iter()) {
      *v = w;
    }
    for (v, &w) in self.interior_loop_length_counts_symm.iter_mut().zip(CONTRA_IL_SYMM_LENGTH_FES_AT_LEAST.iter()) {
      *v = w;
    }
    for (v, &w) in self.interior_loop_length_counts_asymm.iter_mut().zip(CONTRA_IL_ASYMM_LENGTH_FES_AT_LEAST.iter()) {
      *v = w;
    }
    let len = self.stack_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          for l in 0 .. len {
            if !is_canonical(&(k, l)) {continue;}
            self.stack_count_mat[i][j][k][l] = CONTRA_STACK_FES[i][j][k][l];
          }
        }
      }
    }
    let len = self.terminal_mismatch_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          for l in 0 .. len {
            self.terminal_mismatch_count_mat[i][j][k][l] = CONTRA_TERMINAL_MISMATCH_FES[i][j][k][l];
          }
        }
      }
    }
    let len = self.left_dangle_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          self.left_dangle_count_mat[i][j][k] = CONTRA_LEFT_DANGLE_FES[i][j][k];
        }
      }
    }
    let len = self.right_dangle_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          self.right_dangle_count_mat[i][j][k] = CONTRA_RIGHT_DANGLE_FES[i][j][k];
        }
      }
    }
    let len = self.helix_end_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        self.helix_end_count_mat[i][j] = CONTRA_HELIX_CLOSING_FES[i][j];
      }
    }
    let len = self.base_pair_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        self.base_pair_count_mat[i][j] = CONTRA_BASE_PAIR_FES[i][j];
      }
    }
    let len = self.interior_loop_length_count_mat_explicit.len();
    for i in 0 .. len {
      for j in 0 .. len {
        self.interior_loop_length_count_mat_explicit[i][j] = CONTRA_IL_EXPLICIT_FES[i][j];
      }
    }
    for (v, &w) in self.bulge_loop_0x1_length_counts.iter_mut().zip(CONTRA_BL_0X1_FES.iter()) {
      *v = w;
    }
    let len = self.interior_loop_1x1_length_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        self.interior_loop_1x1_length_count_mat[i][j] = CONTRA_IL_1X1_FES[i][j];
      }
    }
    self.multi_loop_base_count = CONTRA_ML_BASE_FE;
    self.multi_loop_basepairing_count = CONTRA_ML_PAIRED_FE;
    self.multi_loop_accessible_baseunpairing_count = CONTRA_ML_UNPAIRED_FE;
    self.external_loop_accessible_basepairing_count = CONTRA_EL_PAIRED_FE;
    self.external_loop_accessible_baseunpairing_count = CONTRA_EL_UNPAIRED_FE;
    self.accumulate();
  }
}

impl<T: Hash> SsPartFuncMats<T> {
  pub fn new(seq_len: usize) -> SsPartFuncMats<T> {
    let neg_inf_mat = vec![vec![NEG_INFINITY; seq_len]; seq_len];
    SsPartFuncMats {
      part_func_mat: vec![vec![0.; seq_len]; seq_len],
      part_func_mat_4_rightmost_base_pairings: neg_inf_mat.clone(),
      part_func_mat_4_rightmost_base_pairings_on_mls: neg_inf_mat.clone(),
      part_func_mat_4_base_pairings: SparsePartFuncMat::<T>::default(),
      part_func_mat_4_base_pairings_accessible: SparsePartFuncMat::<T>::default(),
      part_func_mat_4_ml: neg_inf_mat.clone(),
      part_func_mat_4_at_least_1_base_pairings_on_mls: neg_inf_mat,
    }
  }
}

impl<T: Unsigned + PrimInt + Hash + One> SsFreeEnergyMats<T> {
  pub fn new() -> SsFreeEnergyMats<T> {
    let free_energy_mat = SparseFreeEnergyMat::<T>::default();
    let free_energy_4d_mat = FreeEnergy4dMat::<T>::default();
    SsFreeEnergyMats {
      hl_fe_mat: free_energy_mat.clone(),
      twoloop_fe_4d_mat: free_energy_4d_mat,
      ml_closing_bp_fe_mat: free_energy_mat.clone(),
      accessible_bp_fe_mat: free_energy_mat,
    }
  }
}

pub fn mccaskill_algo<T>(seq: SeqSlice, uses_contra_model: bool, allows_short_hairpins: bool, struct_feature_score_sets: &StructFeatureCountSets,) -> (SparseProbMat<T>, SsFreeEnergyMats<T>)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let seq_len = seq.len();
  let mut ss_free_energy_mats = SsFreeEnergyMats::<T>::new();
  let ss_part_func_mats = if uses_contra_model {
    get_ss_part_func_mats_contra::<T>(seq, seq_len, &mut ss_free_energy_mats, allows_short_hairpins, struct_feature_score_sets)
  } else {
    get_ss_part_func_mats::<T>(seq, seq_len, &mut ss_free_energy_mats)
  };
  let bpp_mat = if uses_contra_model {
    get_base_pairing_prob_mat_contra::<T>(&ss_part_func_mats, seq_len, &ss_free_energy_mats, allows_short_hairpins, struct_feature_score_sets)
  } else {
    get_base_pairing_prob_mat::<T>(&ss_part_func_mats, seq_len, &ss_free_energy_mats)
  };
  (bpp_mat, ss_free_energy_mats)
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
          if long_k - long_i - 1 > MAX_2_LOOP_LEN {break;}
          for l in range(k + T::one(), j).rev() {
            let long_l = l.to_usize().unwrap();
            if long_j - long_l - 1 + long_k - long_i - 1 > MAX_2_LOOP_LEN {break;}
            let accessible_pp = (k, l);
            let long_accessible_pp = (long_k, long_l);
            match ss_part_func_mats.part_func_mat_4_base_pairings.get(&accessible_pp) {
              Some(&part_func) => {
                let twoloop_free_energy = get_2_loop_fe(seq, &long_pp_closing_loop, &long_accessible_pp);
                ss_free_energy_mats.twoloop_fe_4d_mat.insert((i, j, k, l), twoloop_free_energy);
                let term = part_func + twoloop_free_energy;
                logsumexp(&mut sum, term);
              }, None => {},
            }
          }
        }
        let ml_closing_basepairing_fe = get_ml_closing_basepairing_fe(seq, &long_pp_closing_loop);
        logsumexp(&mut sum, ss_part_func_mats.part_func_mat_4_ml[long_i + 1][long_j - 1] + ml_closing_basepairing_fe);
        let ml_or_el_accessible_basepairing_fe = get_ml_or_el_accessible_basepairing_fe(seq, &long_pp_closing_loop, false);
        if sum > NEG_INFINITY {
          ss_free_energy_mats.ml_closing_bp_fe_mat.insert(pp_closing_loop, ml_closing_basepairing_fe);
          ss_free_energy_mats.accessible_bp_fe_mat.insert(pp_closing_loop, ml_or_el_accessible_basepairing_fe);
          ss_part_func_mats.part_func_mat_4_base_pairings.insert(pp_closing_loop, sum);
          let sum = sum + ml_or_el_accessible_basepairing_fe;
          ss_part_func_mats.part_func_mat_4_base_pairings_accessible.insert(pp_closing_loop, sum);
        }
      }
      sum = NEG_INFINITY;
      for k in range_inclusive(i + T::one(), j) {
        let accessible_pp = (i, k);
        match ss_part_func_mats.part_func_mat_4_base_pairings_accessible.get(&accessible_pp) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
      }
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[long_i][long_j] = sum;
      sum = 0.;
      for k in long_i .. long_j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j];
        let part_func = if long_i == 0 && k == 0 {0.} else {ss_part_func_mats.part_func_mat[long_i][k - 1]};
        logsumexp(&mut sum, part_func + ss_part_func_4_rightmost_base_pairings);
      }
      ss_part_func_mats.part_func_mat[long_i][long_j] = sum;
      sum = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[long_i][long_j] + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
      let mut sum_2 = NEG_INFINITY;
      for k in long_i + 1 .. long_j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j] + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
        logsumexp(&mut sum, ss_part_func_4_rightmost_base_pairings);
        logsumexp(&mut sum_2, ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][k - 1] + ss_part_func_4_rightmost_base_pairings);
      }
      ss_part_func_mats.part_func_mat_4_ml[long_i][long_j] = sum_2;
      logsumexp(&mut sum, sum_2);
      ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][long_j] = sum;
    }
  }
  ss_part_func_mats
}

pub fn get_ss_part_func_mats_contra<T>(seq: SeqSlice, seq_len: usize, ss_free_energy_mats: &mut SsFreeEnergyMats<T>, allows_short_hairpins: bool, struct_feature_score_sets: &StructFeatureCountSets,) -> SsPartFuncMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let mut ss_part_func_mats = SsPartFuncMats::<T>::new(seq_len);
  let short_seq_len = T::from_usize(seq_len).unwrap();
  for sub_seq_len in range_inclusive(T::one(), short_seq_len) {
    for i in range_inclusive(T::zero(), short_seq_len - sub_seq_len) {
      let j = i + sub_seq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let pp_closing_loop = (i, j);
      let long_pp_closing_loop = (long_i, long_j);
      let bp_closing_loop = (seq[long_i], seq[long_j]);
      let mut sum = NEG_INFINITY;
      if is_canonical(&bp_closing_loop) {
        if allows_short_hairpins || (!allows_short_hairpins && long_pp_closing_loop.1 - long_pp_closing_loop.0 + 1 >= MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL) {
          if long_j - long_i - 1 <= CONTRA_MAX_LOOP_LEN {
            let hl_fe = get_hl_fe_contra(seq, &long_pp_closing_loop, struct_feature_score_sets);
            ss_free_energy_mats.hl_fe_mat.insert(pp_closing_loop, hl_fe);
            logsumexp(&mut sum, hl_fe);
          }
          for k in range(i + T::one(), j - T::one()) {
            let long_k = k.to_usize().unwrap();
            if long_k - long_i - 1 > CONTRA_MAX_LOOP_LEN {break;}
            for l in range(k + T::one(), j).rev() {
              let long_l = l.to_usize().unwrap();
              if long_j - long_l - 1 + long_k - long_i - 1 > CONTRA_MAX_LOOP_LEN {break;}
              let accessible_pp = (k, l);
              let long_accessible_pp = (long_k, long_l);
              match ss_part_func_mats.part_func_mat_4_base_pairings.get(&accessible_pp) {
                Some(&part_func) => {
                  let twoloop_free_energy = get_2_loop_fe_contra(seq, &long_pp_closing_loop, &long_accessible_pp, struct_feature_score_sets);
                  ss_free_energy_mats.twoloop_fe_4d_mat.insert((i, j, k, l), twoloop_free_energy);
                  let term = part_func + twoloop_free_energy;
                  logsumexp(&mut sum, term);
                }, None => {},
              }
            }
          }
          let coefficient = struct_feature_score_sets.multi_loop_base_count + struct_feature_score_sets.multi_loop_basepairing_count + get_contra_junction_fe_multi(seq, &long_pp_closing_loop, seq_len, false, struct_feature_score_sets);
          logsumexp(&mut sum, ss_part_func_mats.part_func_mat_4_ml[long_i + 1][long_j - 1] + coefficient);
          let fe_multi = get_contra_junction_fe_multi(seq, &(long_pp_closing_loop.1, long_pp_closing_loop.0), seq_len, false, struct_feature_score_sets) + struct_feature_score_sets.base_pair_count_mat[bp_closing_loop.0][bp_closing_loop.1];
          if sum > NEG_INFINITY {
            ss_free_energy_mats.ml_closing_bp_fe_mat.insert(pp_closing_loop, coefficient);
            ss_free_energy_mats.accessible_bp_fe_mat.insert(pp_closing_loop, fe_multi);
            ss_part_func_mats.part_func_mat_4_base_pairings.insert(pp_closing_loop, sum);
            let sum = sum + fe_multi;
            ss_part_func_mats.part_func_mat_4_base_pairings_accessible.insert(pp_closing_loop, sum);
          }
        }
      }
      sum = NEG_INFINITY;
      let mut sum_2 = sum;
      for k in range_inclusive(i + T::one(), j) {
        let accessible_pp = (i, k);
        match ss_part_func_mats.part_func_mat_4_base_pairings_accessible.get(&accessible_pp) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + struct_feature_score_sets.external_loop_accessible_basepairing_count + struct_feature_score_sets.external_loop_accessible_baseunpairing_count * (j - k).to_f32().unwrap());
            logsumexp(&mut sum_2, part_func + struct_feature_score_sets.multi_loop_basepairing_count + struct_feature_score_sets.multi_loop_accessible_baseunpairing_count * (j - k).to_f32().unwrap());
          }, None => {},
        }
      }
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[long_i][long_j] = sum;
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[long_i][long_j] = sum_2;
      sum = struct_feature_score_sets.external_loop_accessible_baseunpairing_count * sub_seq_len.to_f32().unwrap();
      for k in long_i .. long_j {
        let ss_part_func_4_rightmost_base_pairings = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings[k][long_j];
        let part_func = if long_i == 0 && k == 0 {0.} else {ss_part_func_mats.part_func_mat[long_i][k - 1]};
        logsumexp(&mut sum, part_func + ss_part_func_4_rightmost_base_pairings);
      }
      ss_part_func_mats.part_func_mat[long_i][long_j] = sum;
      sum = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[long_i][long_j];
      sum_2 = NEG_INFINITY;
      for k in long_i + 1 .. long_j {
        let ss_part_func_4_rightmost_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[k][long_j];
        logsumexp(&mut sum, ss_part_func_4_rightmost_base_pairings_on_mls + struct_feature_score_sets.multi_loop_accessible_baseunpairing_count * (k - long_i) as FreeEnergy);
        logsumexp(&mut sum_2, ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][k - 1] + ss_part_func_4_rightmost_base_pairings_on_mls);
      }
      ss_part_func_mats.part_func_mat_4_ml[long_i][long_j] = sum_2;
      logsumexp(&mut sum, sum_2);
      ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][long_j] = sum;
    }
  }
  ss_part_func_mats
}

fn get_base_pairing_prob_mat<T>(ss_part_func_mats: &SsPartFuncMats<T>, seq_len: usize, ss_free_energy_mats: &SsFreeEnergyMats<T>) -> SparseProbMat<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let ss_part_func = ss_part_func_mats.part_func_mat[0][seq_len - 1];
  let mut bpp_mat = SparseProbMat::<T>::default();
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
        match ss_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
          Some(&part_func) => {
            let bpp = bpp_mat[&pp_closing_loop];
            let ml_closing_basepairing_fe = ss_free_energy_mats.ml_closing_bp_fe_mat[&pp_closing_loop];
            let coefficient = bpp + ml_closing_basepairing_fe - part_func;
            logsumexp(&mut sum_1, coefficient + ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_j + 1][long_k - 1]);
            logsumexp(&mut sum_2, coefficient);
          }, None => {},
        }
      }
      prob_mat_4_mls_1[long_i][long_j] = sum_1;
      prob_mat_4_mls_2[long_i][long_j] = sum_2;
      let accessible_pp = (i, j);
      match ss_part_func_mats.part_func_mat_4_base_pairings.get(&accessible_pp) {
        Some(&part_func) => {
          let ss_part_func_4_base_pairing_accessible = ss_part_func_mats.part_func_mat_4_base_pairings_accessible[&accessible_pp];
          let part_func_pair = (
            if accessible_pp.0 < T::one() {0.} else {ss_part_func_mats.part_func_mat[0][long_i - 1]},
            if accessible_pp.1 > short_seq_len - T::from_usize(2).unwrap() {0.} else {ss_part_func_mats.part_func_mat[long_j + 1][seq_len - 1]},
          );
          let mut sum = part_func_pair.0 + ss_part_func_4_base_pairing_accessible + part_func_pair.1 - ss_part_func;
          for k in range(T::zero(), i).rev() {
            let long_k = k.to_usize().unwrap();
            if long_i - long_k - 1 > MAX_2_LOOP_LEN {break;}
            for l in range(j + T::one(), short_seq_len) {
              let long_l = l.to_usize().unwrap();
              if long_l - long_j - 1 + long_i - long_k - 1 > MAX_2_LOOP_LEN {break;}
              let pp_closing_loop = (k, l);
              match ss_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
                Some(&part_func_2) => {
                  logsumexp(&mut sum, bpp_mat[&pp_closing_loop] + part_func - part_func_2 + ss_free_energy_mats.twoloop_fe_4d_mat[&(k, l, i, j)]);
                }, None => {},
              }
            }
          }
          let coefficient = ss_part_func_4_base_pairing_accessible + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE;
          for k in 0 .. long_i {
            let ss_part_func_4_at_least_1_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_i - 1];
            logsumexp(&mut sum, coefficient + prob_mat_4_mls_2[k][long_j] + ss_part_func_4_at_least_1_base_pairings_on_mls);
            let prob_4_mls = prob_mat_4_mls_1[k][long_j];
            logsumexp(&mut sum, coefficient + prob_4_mls);
            logsumexp(&mut sum, coefficient + prob_4_mls + ss_part_func_4_at_least_1_base_pairings_on_mls);
          }
          debug_assert!(NEG_INFINITY <= sum && sum <= 0.);
          if sum > NEG_INFINITY {
            bpp_mat.insert(accessible_pp, sum);
          }
        }, None => {},
      }
    }
  }
  bpp_mat = bpp_mat.iter().map(|(pos_pair, &bpp)| {(*pos_pair, expf(bpp))}).collect();
  bpp_mat
}

fn get_base_pairing_prob_mat_contra<T>(ss_part_func_mats: &SsPartFuncMats<T>, seq_len: usize, ss_free_energy_mats: &SsFreeEnergyMats<T>, allows_short_hairpins: bool, struct_feature_score_sets: &StructFeatureCountSets,) -> SparseProbMat<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let ss_part_func = ss_part_func_mats.part_func_mat[0][seq_len - 1];
  let mut bpp_mat = SparseProbMat::<T>::default();
  let mut prob_mat_4_mls_1 = vec![vec![NEG_INFINITY; seq_len]; seq_len];
  let mut prob_mat_4_mls_2 = prob_mat_4_mls_1.clone();
  let short_seq_len = T::from_usize(seq_len).unwrap();
  for sub_seq_len in range_inclusive(T::from_usize(if allows_short_hairpins {2} else {MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL}).unwrap(), short_seq_len).rev() {
    for i in range_inclusive(T::zero(), short_seq_len - sub_seq_len) {
      let j = i + sub_seq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let mut sum_1 = NEG_INFINITY;
      let mut sum_2 = sum_1;
      for k in range(j + T::one(), short_seq_len) {
        let long_k = k.to_usize().unwrap();
        let pp_closing_loop = (i, k);
        match ss_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
          Some(&part_func) => {
            let bpp = bpp_mat[&pp_closing_loop];
            let ml_closing_basepairing_fe = ss_free_energy_mats.ml_closing_bp_fe_mat[&pp_closing_loop];
            let coefficient = bpp + ml_closing_basepairing_fe - part_func;
            logsumexp(&mut sum_1, coefficient + ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_j + 1][long_k - 1]);
            logsumexp(&mut sum_2, coefficient + struct_feature_score_sets.multi_loop_accessible_baseunpairing_count * (k - j - T::one()).to_f32().unwrap());
          }, None => {},
        }
      }
      prob_mat_4_mls_1[long_i][long_j] = sum_1;
      prob_mat_4_mls_2[long_i][long_j] = sum_2;
      let accessible_pp = (i, j);
      match ss_part_func_mats.part_func_mat_4_base_pairings.get(&accessible_pp) {
        Some(&part_func) => {
          let part_func_pair = (
            if accessible_pp.0 < T::one() {0.} else {ss_part_func_mats.part_func_mat[0][long_i - 1]},
            if accessible_pp.1 > short_seq_len - T::from_usize(2).unwrap() {0.} else {ss_part_func_mats.part_func_mat[long_j + 1][seq_len - 1]},
          );
          let mut sum = part_func_pair.0 + part_func_pair.1 + ss_part_func_mats.part_func_mat_4_base_pairings_accessible[&accessible_pp] + struct_feature_score_sets.external_loop_accessible_basepairing_count - ss_part_func;
          for k in range(T::zero(), i).rev() {
            let long_k = k.to_usize().unwrap();
            if long_i - long_k - 1 > CONTRA_MAX_LOOP_LEN {break;}
            for l in range(j + T::one(), short_seq_len) {
              let long_l = l.to_usize().unwrap();
              if long_l - long_j - 1 + long_i - long_k - 1 > CONTRA_MAX_LOOP_LEN {break;}
              let pp_closing_loop = (k, l);
              match ss_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
                Some(&part_func_2) => {
                  logsumexp(&mut sum, bpp_mat[&pp_closing_loop] + part_func - part_func_2 + ss_free_energy_mats.twoloop_fe_4d_mat[&(k, l, i, j)]);
                }, None => {},
              }
            }
          }
          let coefficient = ss_part_func_mats.part_func_mat_4_base_pairings_accessible[&accessible_pp] + struct_feature_score_sets.multi_loop_basepairing_count;
          for k in 0 .. long_i {
            let ss_part_func_4_at_least_1_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_i - 1];
            logsumexp(&mut sum, coefficient + prob_mat_4_mls_2[k][long_j] + ss_part_func_4_at_least_1_base_pairings_on_mls);
            let prob_4_mls = prob_mat_4_mls_1[k][long_j];
            logsumexp(&mut sum, coefficient + prob_4_mls + struct_feature_score_sets.multi_loop_accessible_baseunpairing_count * (long_i - k - 1) as FreeEnergy);
            logsumexp(&mut sum, coefficient + prob_4_mls + ss_part_func_4_at_least_1_base_pairings_on_mls);
          }
          debug_assert!(NEG_INFINITY <= sum && sum <= 0.);
          if sum > NEG_INFINITY {
            bpp_mat.insert(accessible_pp, sum);
          }
        }, None => {},
      }
    }
  }
  bpp_mat = bpp_mat.iter().map(|(pos_pair, &bpp)| {(*pos_pair, expf(bpp))}).collect();
  bpp_mat
}
