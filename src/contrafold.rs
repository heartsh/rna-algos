// pub mod trained_feature_score_sets;

use utils::*;
use mccaskill_algo::*;

pub type HairpinLoopLengthCounts =
  [FeatureCount; CONSPROB_MAX_HAIRPIN_LOOP_LEN + 1];
pub type BulgeLoopLengthCounts = [FeatureCount; CONSPROB_MAX_TWOLOOP_LEN];
pub type InteriorLoopLengthCounts = [FeatureCount; CONSPROB_MAX_TWOLOOP_LEN - 1];
pub type InteriorLoopLengthCountsSymm = [FeatureCount; CONSPROB_MAX_INTERIOR_LOOP_LEN_SYMM];
pub type InteriorLoopLengthCountsAsymm = [FeatureCount; CONSPROB_MAX_INTERIOR_LOOP_LEN_ASYMM];
pub type DangleCount3dMat = [[[FeatureCount; NUM_OF_BASES]; NUM_OF_BASES]; NUM_OF_BASES];
pub type BasePairCountMat = HelixEndCountMat;
pub type InteriorLoopLengthCountMatExplicit = [[FeatureCount; CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT]; CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT];
pub type BulgeLoop0x1LengthCounts = [FeatureCount; NUM_OF_BASES];
pub type InteriorLoop1x1LengthCountMat = [[FeatureCount; NUM_OF_BASES]; NUM_OF_BASES];
pub type InteriorLoopLengthCountMat =
  [[FeatureCount; CONSPROB_MAX_TWOLOOP_LEN - 1]; CONSPROB_MAX_TWOLOOP_LEN - 1];
pub type TerminalMismatchCount3dMat = [[[FeatureCount; NUM_OF_BASES]; NUM_OF_BASES]; NUM_OF_BASES];
pub type TerminalMismatchCount4dMat = [TerminalMismatchCount3dMat; NUM_OF_BASES];
pub type StackCountMat = TerminalMismatchCount4dMat;
pub type HelixEndCountMat = [[FeatureCount; NUM_OF_BASES]; NUM_OF_BASES];
pub type TrainData = Vec<TrainDatum>;

#[derive(Clone)]
pub struct TrainDatum {
  pub seq: Seq,
  pub observed_feature_count_sets: FeatureCountSets,
  pub expected_feature_count_sets: FeatureCountSets,
  pub part_func: Prob,
}

#[derive(Clone, Debug)]
pub struct FeatureCountSets {
  pub hairpin_loop_length_counts: HairpinLoopLengthCounts,
  pub bulge_loop_length_counts: BulgeLoopLengthCounts,
  pub interior_loop_length_counts: InteriorLoopLengthCounts,
  pub interior_loop_length_counts_symm: InteriorLoopLengthCountsSymm,
  pub interior_loop_length_counts_asymm: InteriorLoopLengthCountsAsymm,
  pub stack_count_mat: StackCountMat,
  pub terminal_mismatch_count_mat: TerminalMismatchCount4dMat,
  pub left_dangle_count_mat: DangleCount3dMat,
  pub right_dangle_count_mat: DangleCount3dMat,
  pub helix_end_count_mat: HelixEndCountMat,
  pub base_pair_count_mat: BasePairCountMat,
  pub interior_loop_length_count_mat_explicit: InteriorLoopLengthCountMatExplicit,
  pub bulge_loop_0x1_length_counts: BulgeLoop0x1LengthCounts,
  pub interior_loop_1x1_length_count_mat: InteriorLoop1x1LengthCountMat,
  pub multi_loop_base_count: FeatureCount,
  pub multi_loop_basepairing_count: FeatureCount,
  pub multi_loop_accessible_baseunpairing_count: FeatureCount,
  pub external_loop_accessible_basepairing_count: FeatureCount,
  pub external_loop_accessible_baseunpairing_count: FeatureCount,
  pub hairpin_loop_length_counts_cumulative: HairpinLoopLengthCounts,
  pub bulge_loop_length_counts_cumulative: BulgeLoopLengthCounts,
  pub interior_loop_length_counts_cumulative: InteriorLoopLengthCounts,
  pub interior_loop_length_counts_symm_cumulative: InteriorLoopLengthCountsSymm,
  pub interior_loop_length_counts_asymm_cumulative: InteriorLoopLengthCountsAsymm,
}

impl FeatureCountSets {
  pub fn new(init_val: FeatureCount) -> FeatureCountSets {
    // let init_vals = [init_val; NUM_OF_BASES];
    let twod_mat = [[init_val; NUM_OF_BASES]; NUM_OF_BASES];
    let threed_mat = [[[init_val; NUM_OF_BASES]; NUM_OF_BASES]; NUM_OF_BASES];
    let fourd_mat = [[[[init_val; NUM_OF_BASES]; NUM_OF_BASES]; NUM_OF_BASES]; NUM_OF_BASES];
    FeatureCountSets {
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
      hairpin_loop_length_counts_cumulative: [init_val; CONSPROB_MAX_HAIRPIN_LOOP_LEN + 1],
      bulge_loop_length_counts_cumulative: [init_val; CONSPROB_MAX_TWOLOOP_LEN],
      interior_loop_length_counts_cumulative: [init_val; CONSPROB_MAX_TWOLOOP_LEN - 1],
      interior_loop_length_counts_symm_cumulative: [init_val; CONSPROB_MAX_INTERIOR_LOOP_LEN_SYMM],
      interior_loop_length_counts_asymm_cumulative: [init_val; CONSPROB_MAX_INTERIOR_LOOP_LEN_ASYMM],
    }
  }

  pub fn get_len(&self) -> usize {
    self.hairpin_loop_length_counts.len()
      + self.bulge_loop_length_counts.len()
      + self.interior_loop_length_counts.len()
      + self.interior_loop_length_counts_symm.len()
      + self.interior_loop_length_counts_asymm.len()
      + self.stack_count_mat.len().pow(4)
      + self.terminal_mismatch_count_mat.len().pow(4)
      + self.left_dangle_count_mat.len().pow(3)
      + self.right_dangle_count_mat.len().pow(3)
      + self.helix_end_count_mat.len().pow(2)
      + self.base_pair_count_mat.len().pow(2)
      + self.interior_loop_length_count_mat_explicit.len().pow(2)
      + self.bulge_loop_0x1_length_counts.len()
      + self.interior_loop_1x1_length_count_mat.len().pow(2)
      + 1
      + 1
      + 1
      + 1
      + 1
  }

  pub fn update_regularizers(&self, regularizers: &mut Regularizers) {
    let mut regularizers_tmp = vec![0.; regularizers.len()];
    let mut offset = 0;
    let len = self.hairpin_loop_length_counts.len();
    let group_size = len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      let count = self.hairpin_loop_length_counts[i];
      squared_sum += count * count;
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      regularizers_tmp[offset + i] = regularizer;
    }
    offset += group_size;
    let len = self.bulge_loop_length_counts.len();
    let group_size = len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      let count = self.bulge_loop_length_counts[i];
      squared_sum += count * count;
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      regularizers_tmp[offset + i] = regularizer;
    }
    offset += group_size;
    let len = self.interior_loop_length_counts.len();
    let group_size = len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      let count = self.interior_loop_length_counts[i];
      squared_sum += count * count;
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      regularizers_tmp[offset + i] = regularizer;
    }
    offset += group_size;
    let len = self.interior_loop_length_counts_symm.len();
    let group_size = len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      let count = self.interior_loop_length_counts_symm[i];
      squared_sum += count * count;
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      regularizers_tmp[offset + i] = regularizer;
    }
    offset += group_size;
    let len = self.interior_loop_length_counts_asymm.len();
    let group_size = len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      let count = self.interior_loop_length_counts_asymm[i];
      squared_sum += count * count;
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      regularizers_tmp[offset + i] = regularizer;
    }
    offset += group_size;
    let len = self.stack_count_mat.len();
    let group_size = len.pow(4);
    let effective_group_size = NUM_OF_BASEPAIRINGS * NUM_OF_BASEPAIRINGS;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          for l in 0 .. len {
            if !is_canonical(&(k, l)) {continue;}
            let count = self.stack_count_mat[i][j][k][l];
            squared_sum += count * count;
          }
        }
      }
    }
    let regularizer = get_regularizer(effective_group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          for l in 0 .. len {
            if !is_canonical(&(k, l)) {continue;}
            regularizers_tmp[offset + i * len.pow(3) + j * len.pow(2) + k * len + l] = regularizer;
          }
        }
      }
    }
    offset += group_size;
    let len = self.terminal_mismatch_count_mat.len();
    let group_size = len.pow(4);
    let effective_group_size = NUM_OF_BASEPAIRINGS * len * len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          for l in 0 .. len {
            let count = self.terminal_mismatch_count_mat[i][j][k][l];
            squared_sum += count * count;
          }
        }
      }
    }
    let regularizer = get_regularizer(effective_group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          for l in 0 .. len {
            regularizers_tmp[offset + i * len.pow(3) + j * len.pow(2) + k * len + l] = regularizer;
          }
        }
      }
    }
    offset += group_size;
    let len = self.left_dangle_count_mat.len();
    let group_size = len.pow(3);
    let effective_group_size = NUM_OF_BASEPAIRINGS * len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          let count = self.left_dangle_count_mat[i][j][k];
          squared_sum += count * count;
        }
      }
    }
    let regularizer = get_regularizer(effective_group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          regularizers_tmp[offset + i * len.pow(2) + j * len + k] = regularizer;
        }
      }
    }
    offset += group_size;
    let len = self.right_dangle_count_mat.len();
    let group_size = len.pow(3);
    let effective_group_size = NUM_OF_BASEPAIRINGS * len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          let count = self.right_dangle_count_mat[i][j][k];
          squared_sum += count * count;
        }
      }
    }
    let regularizer = get_regularizer(effective_group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          regularizers_tmp[offset + i * len.pow(2) + j * len + k] = regularizer;
        }
      }
    }
    offset += group_size;
    let len = self.helix_end_count_mat.len();
    let group_size = len.pow(2);
    let effective_group_size = NUM_OF_BASEPAIRINGS;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        let count = self.helix_end_count_mat[i][j];
        squared_sum += count * count;
      }
    }
    let regularizer = get_regularizer(effective_group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        regularizers_tmp[offset + i * len + j] = regularizer;
      }
    }
    offset += group_size;
    let len = self.base_pair_count_mat.len();
    let group_size = len.pow(2);
    let effective_group_size = NUM_OF_BASEPAIRINGS;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        let count = self.base_pair_count_mat[i][j];
        squared_sum += count * count;
      }
    }
    let regularizer = get_regularizer(effective_group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        regularizers_tmp[offset + i * len + j] = regularizer;
      }
    }
    offset += group_size;
    let len = self.interior_loop_length_count_mat_explicit.len();
    let group_size = len.pow(2);
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        let count = self.interior_loop_length_count_mat_explicit[i][j];
        squared_sum += count * count;
      }
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        regularizers_tmp[offset + i * len + j] = regularizer;
      }
    }
    offset += group_size;
    let len = self.bulge_loop_0x1_length_counts.len();
    let group_size = len;
    let mut squared_sum = 0.;
    for i in 0 .. len {
      let count = self.bulge_loop_0x1_length_counts[i];
      squared_sum += count * count;
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      regularizers_tmp[offset + i] = regularizer;
    }
    offset += group_size;
    let len = self.interior_loop_1x1_length_count_mat.len();
    let group_size = len.pow(2);
    let mut squared_sum = 0.;
    for i in 0 .. len {
      for j in 0 .. len {
        let count = self.interior_loop_1x1_length_count_mat[i][j];
        squared_sum += count * count;
      }
    }
    let regularizer = get_regularizer(group_size, squared_sum);
    for i in 0 .. len {
      for j in 0 .. len {
        regularizers_tmp[offset + i * len + j] = regularizer;
      }
    }
    offset += group_size;
    let regularizer = get_regularizer(1, self.multi_loop_base_count * self.multi_loop_base_count);
    regularizers_tmp[offset] = regularizer;
    offset += 1;
    let regularizer = get_regularizer(1, self.multi_loop_basepairing_count * self.multi_loop_basepairing_count);
    regularizers_tmp[offset] = regularizer;
    offset += 1;
    let regularizer = get_regularizer(1, self.multi_loop_accessible_baseunpairing_count * self.multi_loop_accessible_baseunpairing_count);
    regularizers_tmp[offset] = regularizer;
    offset += 1;
    let regularizer = get_regularizer(1, self.external_loop_accessible_basepairing_count * self.external_loop_accessible_basepairing_count);
    regularizers_tmp[offset] = regularizer;
    offset += 1;
    let regularizer = get_regularizer(1, self.external_loop_accessible_baseunpairing_count * self.external_loop_accessible_baseunpairing_count);
    regularizers_tmp[offset] = regularizer;
    *regularizers = Array1::from(regularizers_tmp);
  }

  pub fn update(&mut self, train_data: &[TrainDatum], regularizers: &mut Regularizers)
  {
    let f = |_: &BfgsFeatureCounts| {
      self.get_cost(&train_data[..], regularizers) as BfgsFeatureCount
    };
    let g = |_: &BfgsFeatureCounts| {
      convert_feature_counts_2_bfgs_feature_counts(&self.get_grad(train_data, regularizers))
    };
    match bfgs(convert_feature_counts_2_bfgs_feature_counts(&convert_struct_2_vec(self, false)), f, g) {
      Ok(solution) => {
        *self = convert_vec_2_struct(&convert_bfgs_feature_counts_2_feature_counts(&solution), false);
      }, Err(_) => {
        println!("BFGS failed");
      },
    };
    self.update_regularizers(regularizers);
    self.accumulate();
    let vec = convert_struct_2_vec(self, false);
    println!("norm: {}", (&vec * &vec).sum());
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

  pub fn get_grad(&self, train_data: &[TrainDatum], regularizers: &Regularizers) -> FeatureCounts
  {
    let feature_scores = convert_struct_2_vec(self, false);
    let mut grad = FeatureCountSets::new(0.);
    for train_datum in train_data {
      let ref obs = train_datum.observed_feature_count_sets;
      let ref expect = train_datum.expected_feature_count_sets;
      for i in 0 .. obs.hairpin_loop_length_counts.len() {
        for j in 0 .. i + 1 {
          grad.hairpin_loop_length_counts[j] -= obs.hairpin_loop_length_counts[i] - expect.hairpin_loop_length_counts[i];
        }
      }
      for i in 0 .. obs.bulge_loop_length_counts.len() {
        for j in 0 .. i + 1 {
          grad.bulge_loop_length_counts[j] -= obs.bulge_loop_length_counts[i] - expect.bulge_loop_length_counts[i];
        }
      }
      let len = obs.interior_loop_length_counts.len();
      for i in 0 .. len {
        let obs_count = obs.interior_loop_length_counts[i];
        let expect_count = expect.interior_loop_length_counts[i];
        for j in 0 .. i + 1 {
          grad.interior_loop_length_counts[j] -= obs_count - expect_count;
        }
      }
      let len = obs.interior_loop_length_counts_symm.len();
      for i in 0 .. len {
        let obs_count = obs.interior_loop_length_counts_symm[i];
        let expect_count = expect.interior_loop_length_counts_symm[i];
        for j in 0 .. i + 1 {
          grad.interior_loop_length_counts_symm[j] -= obs_count - expect_count;
        }
      }
      let len = obs.interior_loop_length_counts_asymm.len();
      for i in 0 .. len {
        let obs_count = obs.interior_loop_length_counts_asymm[i];
        let expect_count = expect.interior_loop_length_counts_asymm[i];
        for j in 0 .. i + 1 {
          grad.interior_loop_length_counts_asymm[j] -= obs_count - expect_count;
        }
      }
      for i in 0 .. NUM_OF_BASES {
        for j in 0 .. NUM_OF_BASES {
          if !is_canonical(&(i, j)) {continue;}
          for k in 0 .. NUM_OF_BASES {
            for l in 0 .. NUM_OF_BASES {
              if !is_canonical(&(k, l)) {continue;}
              let dict_min_stack = get_dict_min_stack(&(i, j), &(k, l));
              let obs_count = obs.stack_count_mat[dict_min_stack.0.0][dict_min_stack.0.1][dict_min_stack.1.0][dict_min_stack.1.1];
              let expect_count = expect.stack_count_mat[dict_min_stack.0.0][dict_min_stack.0.1][dict_min_stack.1.0][dict_min_stack.1.1];
              grad.stack_count_mat[i][j][k][l] -= obs_count - expect_count;
            }
          }
        }
      }
      for i in 0 .. NUM_OF_BASES {
        for j in 0 .. NUM_OF_BASES {
          if !is_canonical(&(i, j)) {continue;}
          for k in 0 .. NUM_OF_BASES {
            for l in 0 .. NUM_OF_BASES {
              let obs_count = obs.terminal_mismatch_count_mat[i][j][k][l];
              let expect_count = expect.terminal_mismatch_count_mat[i][j][k][l];
              grad.terminal_mismatch_count_mat[i][j][k][l] -= obs_count - expect_count;
            }
          }
        }
      }
      for i in 0 .. NUM_OF_BASES {
        for j in 0 .. NUM_OF_BASES {
          if !is_canonical(&(i, j)) {continue;}
          for k in 0 .. NUM_OF_BASES {
            let obs_count = obs.left_dangle_count_mat[i][j][k];
            let expect_count = expect.left_dangle_count_mat[i][j][k];
            grad.left_dangle_count_mat[i][j][k] -= obs_count - expect_count;
          }
        }
      }
      for i in 0 .. NUM_OF_BASES {
        for j in 0 .. NUM_OF_BASES {
          if !is_canonical(&(i, j)) {continue;}
          for k in 0 .. NUM_OF_BASES {
            let obs_count = obs.right_dangle_count_mat[i][j][k];
            let expect_count = expect.right_dangle_count_mat[i][j][k];
            grad.right_dangle_count_mat[i][j][k] -= obs_count - expect_count;
          }
        }
      }
      for i in 0 .. NUM_OF_BASES {
        for j in 0 .. NUM_OF_BASES {
          if !is_canonical(&(i, j)) {continue;}
          let obs_count = obs.helix_end_count_mat[i][j];
          let expect_count = expect.helix_end_count_mat[i][j];
          grad.helix_end_count_mat[i][j] -= obs_count - expect_count;
        }
      }
      for i in 0 .. NUM_OF_BASES {
        for j in 0 .. NUM_OF_BASES {
          if !is_canonical(&(i, j)) {continue;}
          let dict_min_base_pair = get_dict_min_base_pair(&(i, j));
          let obs_count = obs.base_pair_count_mat[dict_min_base_pair.0][dict_min_base_pair.1];
          let expect_count = expect.base_pair_count_mat[dict_min_base_pair.0][dict_min_base_pair.1];
          grad.base_pair_count_mat[i][j] -= obs_count - expect_count;
        }
      }
      let len = obs.interior_loop_length_count_mat_explicit.len();
      for i in 0 .. len {
        for j in 0 .. len {
          let dict_min_loop_len_pair = get_dict_min_loop_len_pair(&(i, j));
          let obs_count = obs.interior_loop_length_count_mat_explicit[dict_min_loop_len_pair.0][dict_min_loop_len_pair.1];
          let expect_count = expect.interior_loop_length_count_mat_explicit[dict_min_loop_len_pair.0][dict_min_loop_len_pair.1];
          grad.interior_loop_length_count_mat_explicit[i][j] -= obs_count - expect_count;
        }
      }
      for i in 0 .. NUM_OF_BASES {
        let obs_count = obs.bulge_loop_0x1_length_counts[i];
        let expect_count = expect.bulge_loop_0x1_length_counts[i];
        grad.bulge_loop_0x1_length_counts[i] -= obs_count - expect_count;
      }
      for i in 0 .. NUM_OF_BASES {
        for j in 0 .. NUM_OF_BASES {
          let dict_min_nuc_pair = get_dict_min_nuc_pair(&(i, j));
          let obs_count = obs.interior_loop_1x1_length_count_mat[dict_min_nuc_pair.0][dict_min_nuc_pair.1];
          let expect_count = expect.interior_loop_1x1_length_count_mat[dict_min_nuc_pair.0][dict_min_nuc_pair.1];
          grad.interior_loop_1x1_length_count_mat[i][j] -= obs_count - expect_count;
        }
      }
      let obs_count = obs.multi_loop_base_count;
      let expect_count = expect.multi_loop_base_count;
      grad.multi_loop_base_count -= obs_count - expect_count;
      let obs_count = obs.multi_loop_basepairing_count;
      let expect_count = expect.multi_loop_basepairing_count;
      grad.multi_loop_basepairing_count -= obs_count - expect_count;
      let obs_count = obs.multi_loop_accessible_baseunpairing_count;
      let expect_count = expect.multi_loop_accessible_baseunpairing_count;
      grad.multi_loop_accessible_baseunpairing_count -= obs_count - expect_count;
      let obs_count = obs.external_loop_accessible_basepairing_count;
      let expect_count = expect.external_loop_accessible_basepairing_count;
      grad.external_loop_accessible_basepairing_count -= obs_count - expect_count;
      let obs_count = obs.external_loop_accessible_baseunpairing_count;
      let expect_count = expect.external_loop_accessible_baseunpairing_count;
      grad.external_loop_accessible_baseunpairing_count -= obs_count - expect_count;
    }
    convert_struct_2_vec(&grad, false) + regularizers.clone() * feature_scores
  }

  pub fn get_cost(&self, train_data: &[TrainDatum], regularizers: &Regularizers) -> FeatureCount {
    let mut log_likelihood = 0.;
    let feature_scores_cumulative = convert_struct_2_vec(self, true);
    for train_datum in train_data {
      let ref obs = train_datum.observed_feature_count_sets;
      log_likelihood += feature_scores_cumulative.dot(&convert_struct_2_vec(obs, false));
      log_likelihood -= train_datum.part_func;
    }
    let feature_scores = convert_struct_2_vec(self, false);
    let product = regularizers.clone() * feature_scores.clone();
    - log_likelihood + product.dot(&feature_scores) / 2.
  }

  pub fn rand_init(&mut self) {
    let len = self.get_len();
    let std_deviation = 1. / (len as FeatureCount).sqrt();
    let normal = Normal::new(0., std_deviation).unwrap();
    let mut thread_rng = thread_rng();
    for v in self.hairpin_loop_length_counts.iter_mut() {
      *v = normal.sample(&mut thread_rng);
    }
    for v in self.bulge_loop_length_counts.iter_mut() {
      *v = normal.sample(&mut thread_rng);
    }
    for v in self.interior_loop_length_counts.iter_mut() {
      *v = normal.sample(&mut thread_rng);
    }
    for v in self.interior_loop_length_counts_symm.iter_mut() {
      *v = normal.sample(&mut thread_rng);
    }
    for v in self.interior_loop_length_counts_asymm.iter_mut() {
      *v = normal.sample(&mut thread_rng);
    }
    let len = self.stack_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          for l in 0 .. len {
            if !is_canonical(&(k, l)) {continue;}
            let v = normal.sample(&mut thread_rng);
            if self.stack_count_mat[i][j][k][l] == 0. {
              self.stack_count_mat[i][j][k][l] = v;
            }
            if self.stack_count_mat[l][k][j][i] == 0. {
              self.stack_count_mat[l][k][j][i] = v;
            }
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
            let v = normal.sample(&mut thread_rng);
            self.terminal_mismatch_count_mat[i][j][k][l] = v;
          }
        }
      }
    }
    let len = self.left_dangle_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          let v = normal.sample(&mut thread_rng);
          self.left_dangle_count_mat[i][j][k] = v;
        }
      }
    }
    let len = self.right_dangle_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        for k in 0 .. len {
          let v = normal.sample(&mut thread_rng);
          self.right_dangle_count_mat[i][j][k] = v;
        }
      }
    }
    let len = self.helix_end_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        let v = normal.sample(&mut thread_rng);
        self.helix_end_count_mat[i][j] = v;
      }
    }
    let len = self.base_pair_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        if !is_canonical(&(i, j)) {continue;}
        let v = normal.sample(&mut thread_rng);
        if self.base_pair_count_mat[i][j] == 0. {
          self.base_pair_count_mat[i][j] = v;
        }
        if self.base_pair_count_mat[j][i] == 0. {
          self.base_pair_count_mat[j][i] = v;
        }
      }
    }
    let len = self.interior_loop_length_count_mat_explicit.len();
    for i in 0 .. len {
      for j in 0 .. len {
        let v = normal.sample(&mut thread_rng);
        if self.interior_loop_length_count_mat_explicit[i][j] == 0. {
          self.interior_loop_length_count_mat_explicit[i][j] = v;
        }
        if self.interior_loop_length_count_mat_explicit[j][i] == 0. {
          self.interior_loop_length_count_mat_explicit[j][i] = v;
        }
      }
    }
    for v in &mut self.bulge_loop_0x1_length_counts {
      *v = normal.sample(&mut thread_rng);
    }
    let len = self.interior_loop_1x1_length_count_mat.len();
    for i in 0 .. len {
      for j in 0 .. len {
        let v = normal.sample(&mut thread_rng);
        if self.interior_loop_1x1_length_count_mat[i][j] == 0. {
          self.interior_loop_1x1_length_count_mat[i][j] = v;
        }
        if self.interior_loop_1x1_length_count_mat[j][i] == 0. {
          self.interior_loop_1x1_length_count_mat[j][i] = v;
        }
      }
    }
    self.multi_loop_base_count = normal.sample(&mut thread_rng);
    self.multi_loop_basepairing_count = normal.sample(&mut thread_rng);
    self.multi_loop_accessible_baseunpairing_count = normal.sample(&mut thread_rng);
    self.external_loop_accessible_basepairing_count = normal.sample(&mut thread_rng);
    self.external_loop_accessible_baseunpairing_count = normal.sample(&mut thread_rng);
    self.accumulate();
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

impl TrainDatum {
  pub fn origin() -> TrainDatum {
    TrainDatum {
      seq: Seq::new(),
      observed_feature_count_sets: FeatureCountSets::new(0.),
      expected_feature_count_sets: FeatureCountSets::new(NEG_INFINITY),
      part_func: NEG_INFINITY,
    }
  }

  pub fn new(input_file_path: &Path) -> TrainDatum {
    let fasta_file_reader = Reader::from_file(Path::new(input_file_path)).unwrap();
    let fasta_records: Vec<Record> = fasta_file_reader.records().map(|rec| {rec.unwrap()}).collect();
    let mut seq = convert(fasta_records[0].seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    let second_struct = fasta_records[1].seq();
    // second_struct.insert(0, UNPAIRING_BASE);
    // second_struct.push(UNPAIRING_BASE);
    let second_struct = [&[UNPAIRING_BASE], second_struct, &[UNPAIRING_BASE]].concat();
    let mut train_datum = TrainDatum {
      seq: seq,
      observed_feature_count_sets: FeatureCountSets::new(0.),
      expected_feature_count_sets: FeatureCountSets::new(NEG_INFINITY),
      part_func: NEG_INFINITY,
    };
    train_datum.convert(&second_struct);
    train_datum
  }

  pub fn convert(&mut self, dot_bracket_notation: TextSlice) {
    let seq_len = self.seq.len();
    let mut stack = Vec::new();
    let mut second_struct = HashSet::<(usize, usize)>::default();
    // for i in 0 .. seq_len {
    for i in 1 .. seq_len - 1 {
      let notation_char = dot_bracket_notation[i];
      if notation_char == BASE_PAIRING_LEFT_BASE {
        stack.push(i);
      } else if notation_char == BASE_PAIRING_RIGHT_BASE {
        let pos = stack.pop().unwrap();
        let base_pair = (self.seq[pos], self.seq[i]);
        if !is_canonical(&base_pair) {continue;}
        second_struct.insert((pos, i));
        let dict_min_base_pair = get_dict_min_base_pair(&base_pair);
        self.observed_feature_count_sets.base_pair_count_mat[dict_min_base_pair.0][dict_min_base_pair.1] += 1.;
      }
    }
    let mut loop_struct = HashMap::<(usize, usize), Vec<(usize, usize)>>::default();
    let mut stored_pos_pairs = HashSet::<(usize, usize)>::default();
    for substr_len in 2 .. seq_len - 1 {
      for i in 1 .. seq_len - substr_len {
        let mut is_pos_pair_found = false;
        let j = i + substr_len - 1;
        if second_struct.contains(&(i, j)) {
          is_pos_pair_found = true;
          second_struct.remove(&(i, j));
          stored_pos_pairs.insert((i, j));
        }
        if is_pos_pair_found {
          let mut pos_pairs_in_loop = Vec::new();
          for stored_pos_pair in stored_pos_pairs.iter() {
            if i < stored_pos_pair.0 && stored_pos_pair.1 < j {
              pos_pairs_in_loop.push(*stored_pos_pair);
            }
          }
          for pos_pair in pos_pairs_in_loop.iter() {
            stored_pos_pairs.remove(pos_pair);
          }
          pos_pairs_in_loop.sort();
          loop_struct.insert((i, j), pos_pairs_in_loop);
        }
      }
    }
    for (pos_pair_closing_loop, pos_pairs_in_loop) in loop_struct.iter() {
      let num_of_basepairings_in_loop = pos_pairs_in_loop.len();
      let base_pair = (self.seq[pos_pair_closing_loop.0], self.seq[pos_pair_closing_loop.1]);
      let mismatch_pair = get_mismatch_pair(&self.seq[..], &pos_pair_closing_loop, true);
      if num_of_basepairings_in_loop == 0 {
        self.observed_feature_count_sets.terminal_mismatch_count_mat[base_pair.0][base_pair.1][mismatch_pair.0][mismatch_pair.1] += 1.;
        let hairpin_loop_length = pos_pair_closing_loop.1 - pos_pair_closing_loop.0 - 1;
        if hairpin_loop_length <= CONSPROB_MAX_HAIRPIN_LOOP_LEN {
          self.observed_feature_count_sets.hairpin_loop_length_counts[hairpin_loop_length] += 1.;
        }
        self.observed_feature_count_sets.helix_end_count_mat[base_pair.0][base_pair.1] += 1.;
      } else if num_of_basepairings_in_loop == 1 {
        let ref pos_pair_in_loop = pos_pairs_in_loop[0];
        let base_pair_2 = (self.seq[pos_pair_in_loop.0], self.seq[pos_pair_in_loop.1]);
        let twoloop_length_pair = (pos_pair_in_loop.0 - pos_pair_closing_loop.0 - 1, pos_pair_closing_loop.1 - pos_pair_in_loop.1 - 1);
        let sum = twoloop_length_pair.0 + twoloop_length_pair.1;
        if sum == 0 {
          let dict_min_stack = get_dict_min_stack(&base_pair, &base_pair_2);
          self.observed_feature_count_sets.stack_count_mat[dict_min_stack.0.0][dict_min_stack.0.1][dict_min_stack.1.0][dict_min_stack.1.1] += 1.;
        } else {
          if twoloop_length_pair.0 == 0 || twoloop_length_pair.1 == 0 {
            if sum <= CONSPROB_MAX_TWOLOOP_LEN {
              self.observed_feature_count_sets.bulge_loop_length_counts[sum - 1] += 1.;
              if sum == 1 {
                let mismatch = if twoloop_length_pair.0 == 0 {mismatch_pair.1} else {mismatch_pair.0};
                self.observed_feature_count_sets.bulge_loop_0x1_length_counts[mismatch] += 1.;
              }
            }
          } else {
            if sum <= CONSPROB_MAX_TWOLOOP_LEN {
              self.observed_feature_count_sets.interior_loop_length_counts[sum - 2] += 1.;
              let diff = get_diff(twoloop_length_pair.0, twoloop_length_pair.1);
              if diff == 0 {
                self.observed_feature_count_sets.interior_loop_length_counts_symm[twoloop_length_pair.0 - 1] += 1.;
              } else {
                self.observed_feature_count_sets.interior_loop_length_counts_asymm[diff - 1] += 1.;
              }
              if twoloop_length_pair.0 == 1 && twoloop_length_pair.1 == 1 {
                let dict_min_mismatch_pair = get_dict_min_mismatch_pair(&mismatch_pair);
                self.observed_feature_count_sets.interior_loop_1x1_length_count_mat[dict_min_mismatch_pair.0][dict_min_mismatch_pair.1] += 1.;
              }
              if twoloop_length_pair.0 <= CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT && twoloop_length_pair.1 <= CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT {
                let dict_min_loop_len_pair = get_dict_min_loop_len_pair(&twoloop_length_pair);
                self.observed_feature_count_sets.interior_loop_length_count_mat_explicit[dict_min_loop_len_pair.0 - 1][dict_min_loop_len_pair.1 - 1] += 1.;
              }
            }
          }
          self.observed_feature_count_sets.terminal_mismatch_count_mat[base_pair.0][base_pair.1][mismatch_pair.0][mismatch_pair.1] += 1.;
          let mismatch_pair_2 = get_mismatch_pair(&self.seq[..], pos_pair_in_loop, false);
          self.observed_feature_count_sets.terminal_mismatch_count_mat[base_pair_2.1][base_pair_2.0][mismatch_pair_2.1][mismatch_pair_2.0] += 1.;
          self.observed_feature_count_sets.helix_end_count_mat[base_pair.0][base_pair.1] += 1.;
          self.observed_feature_count_sets.helix_end_count_mat[base_pair_2.1][base_pair_2.0] += 1.;
        }
      } else {
        self.observed_feature_count_sets.left_dangle_count_mat[base_pair.0][base_pair.1][mismatch_pair.0] += 1.;
        self.observed_feature_count_sets.right_dangle_count_mat[base_pair.0][base_pair.1][mismatch_pair.1] += 1.;
        for pos_pair_in_loop in pos_pairs_in_loop.iter() {
          let base_pair_2 = (self.seq[pos_pair_in_loop.0], self.seq[pos_pair_in_loop.1]);
          let mismatch_pair_2 = get_mismatch_pair(&self.seq[..], pos_pair_in_loop, false);
          self.observed_feature_count_sets.left_dangle_count_mat[base_pair_2.1][base_pair_2.0][mismatch_pair_2.1] += 1.;
          self.observed_feature_count_sets.right_dangle_count_mat[base_pair_2.1][base_pair_2.0][mismatch_pair_2.0] += 1.;
          self.observed_feature_count_sets.helix_end_count_mat[base_pair_2.1][base_pair_2.0] += 1.;
        }
        self.observed_feature_count_sets.multi_loop_base_count += 1.;
        self.observed_feature_count_sets.multi_loop_basepairing_count += 1.;
        self.observed_feature_count_sets.multi_loop_basepairing_count += num_of_basepairings_in_loop as Prob;
        let num_of_baseunpairing_nucs_in_multi_loop = get_num_of_multiloop_baseunpairing_nucs(pos_pair_closing_loop, pos_pairs_in_loop);
        self.observed_feature_count_sets.multi_loop_accessible_baseunpairing_count += num_of_baseunpairing_nucs_in_multi_loop as Prob;
        self.observed_feature_count_sets.helix_end_count_mat[base_pair.0][base_pair.1] += 1.;
      }
    }
    self.observed_feature_count_sets.external_loop_accessible_basepairing_count += stored_pos_pairs.len() as Prob;
    let mut stored_pos_pairs_sorted = stored_pos_pairs.iter().map(|x| *x).collect::<Vec<(usize, usize)>>();
    stored_pos_pairs_sorted.sort();
    let num_of_baseunpairing_nucs_in_external_loop = get_num_of_externalloop_baseunpairing_nucs(&stored_pos_pairs_sorted, &self.seq[..]);
    self.observed_feature_count_sets.external_loop_accessible_baseunpairing_count += num_of_baseunpairing_nucs_in_external_loop as Prob;
    for pos_pair_in_loop in stored_pos_pairs.iter() {
      let base_pair = (self.seq[pos_pair_in_loop.0], self.seq[pos_pair_in_loop.1]);
      let mismatch_pair = get_mismatch_pair(&self.seq[..], &pos_pair_in_loop, false);
      if mismatch_pair.1 != PSEUDO_BASE {
        self.observed_feature_count_sets.left_dangle_count_mat[base_pair.1][base_pair.0][mismatch_pair.1] += 1.;
      }
      if mismatch_pair.0 != PSEUDO_BASE {
        self.observed_feature_count_sets.right_dangle_count_mat[base_pair.1][base_pair.0][mismatch_pair.0] += 1.;
      }
      self.observed_feature_count_sets.helix_end_count_mat[base_pair.1][base_pair.0] += 1.;
    }
  }
}

pub const NUM_OF_BASEPAIRINGS: usize = 6;
pub const CONSPROB_MAX_HAIRPIN_LOOP_LEN: usize = 30;
pub const CONSPROB_MAX_TWOLOOP_LEN: usize = CONSPROB_MAX_HAIRPIN_LOOP_LEN;
pub const CONSPROB_MIN_HAIRPIN_LOOP_LEN: usize = 3;
pub const CONSPROB_MIN_HAIRPIN_LOOP_SPAN: usize = CONSPROB_MIN_HAIRPIN_LOOP_LEN + 2;
pub const CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT: usize = 4;
pub const CONSPROB_MAX_INTERIOR_LOOP_LEN_SYMM: usize = CONSPROB_MAX_TWOLOOP_LEN / 2;
pub const CONSPROB_MAX_INTERIOR_LOOP_LEN_ASYMM: usize = CONSPROB_MAX_TWOLOOP_LEN - 2;
pub const TRAINED_FEATURE_SCORE_SETS_FILE_PATH: &'static str = "../src/trained_feature_score_sets_ss.rs";

pub fn contrafold_train<T>(
  thread_pool: &mut Pool,
  train_data: &mut TrainData,
  output_file_path: &Path,
)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send + Display,
{
  let mut feature_score_sets = FeatureCountSets::new(0.);
  feature_score_sets.rand_init();
  // feature_score_sets.transfer();
  let mut old_feature_score_sets = feature_score_sets.clone();
  let mut old_cost = INFINITY;
  let mut costs = Probs::new();
  let mut count = 0;
  let mut regularizers = Regularizers::from(vec![1.; feature_score_sets.get_len()]);
  let num_of_data = train_data.len() as FeatureCount;
  loop {
    let ref ref_2_feature_score_sets = feature_score_sets;
    thread_pool.scoped(|scope| {
      for train_datum in train_data.iter_mut() {
        train_datum.expected_feature_count_sets = FeatureCountSets::new(NEG_INFINITY);
        let seq = &train_datum.seq[..];
        let ref mut expected_feature_count_sets = train_datum.expected_feature_count_sets;
        let ref mut part_func = train_datum.part_func;
        scope.execute(move || {
          *part_func = mccaskill_algo_trained::<T>(
            seq,
            ref_2_feature_score_sets,
            true,
            true,
            expected_feature_count_sets,
          ).1;
        });
      }
    });
    feature_score_sets.update(&train_data, &mut regularizers);
    let cost = feature_score_sets.get_cost(&train_data[..], &regularizers);
    let avg_cost_update_amount = (old_cost - cost) / num_of_data;
    if avg_cost_update_amount < 0. {
      feature_score_sets = old_feature_score_sets.clone();
      break;
    }
    costs.push(cost);
    old_feature_score_sets = feature_score_sets.clone();
    old_cost = cost;
    println!("Epoch {} finished (current cost = {}, average cost update amount = {})", count + 1, cost, avg_cost_update_amount);
    count += 1;
    if avg_cost_update_amount <= LEARNING_TOLERANCE {
      break;
    }
  }
  write_feature_score_sets_trained(&feature_score_sets);
  write_costs(&costs, output_file_path);
}

pub fn mccaskill_algo_trained<T>(seq: SeqSlice, feature_score_sets: &FeatureCountSets, trains_score_params: bool, allows_short_hairpins: bool, expected_feature_count_sets: &mut FeatureCountSets,) -> (SparseProbMat<T>, PartFunc)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let mut ss_free_energy_mats = SsFreeEnergyMats::<T>::new();
  let ss_part_func_mats = get_ss_part_func_mats_trained::<T>(seq, feature_score_sets, &mut ss_free_energy_mats, allows_short_hairpins);
  get_base_pairing_prob_mat_trained::<T>(seq, feature_score_sets, trains_score_params, &ss_part_func_mats, &ss_free_energy_mats, allows_short_hairpins, expected_feature_count_sets)
}

pub fn get_ss_part_func_mats_trained<T>(seq: SeqSlice, feature_score_sets: &FeatureCountSets, ss_free_energy_mats: &mut SsFreeEnergyMats<T>, allows_short_hairpins: bool,) -> SsPartFuncMatsContra<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let seq_len = seq.len();
  let mut ss_part_func_mats = SsPartFuncMatsContra::<T>::new(seq_len);
  let short_seq_len = T::from_usize(seq_len).unwrap();
  for sub_seq_len in range_inclusive(T::one(), short_seq_len - T::from_usize(2).unwrap()) {
    for i in range(T::one(), short_seq_len - sub_seq_len) {
      let j = i + sub_seq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let pp_closing_loop = (i, j);
      let long_pp_closing_loop = (long_i, long_j);
      let bp_closing_loop = (seq[long_i], seq[long_j]);
      let mut sum = NEG_INFINITY;
      if is_canonical(&bp_closing_loop) {
        if !allows_short_hairpins && long_pp_closing_loop.1 - long_pp_closing_loop.0 + 1 < MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL {
          continue;
        }
        if long_j - long_i - 1 <= CONTRA_MAX_LOOP_LEN {
          let hl_fe = get_hl_fe_trained(feature_score_sets, seq, &long_pp_closing_loop);
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
                let twoloop_free_energy = get_2_loop_fe_trained(feature_score_sets, seq, &long_pp_closing_loop, &long_accessible_pp);
                ss_free_energy_mats.twoloop_fe_4d_mat.insert((i, j, k, l), twoloop_free_energy);
                let term = part_func + twoloop_free_energy;
                logsumexp(&mut sum, term);
              }, None => {},
            }
          }
        }
        let coefficient = feature_score_sets.multi_loop_base_count + feature_score_sets.multi_loop_basepairing_count + get_junction_fe_multi_trained(feature_score_sets, seq, &long_pp_closing_loop, seq_len);
        ss_free_energy_mats.ml_closing_bp_fe_mat.insert(pp_closing_loop, coefficient);
        logsumexp(&mut sum, ss_part_func_mats.part_func_mat_4_ml[long_i + 1][long_j - 1] + coefficient);
        let fe_multi = get_junction_fe_multi_trained(feature_score_sets, seq, &(long_pp_closing_loop.1, long_pp_closing_loop.0), seq_len) + feature_score_sets.base_pair_count_mat[bp_closing_loop.0][bp_closing_loop.1];
        ss_free_energy_mats.accessible_bp_fe_mat.insert(pp_closing_loop, fe_multi);
        if sum > NEG_INFINITY {
          ss_part_func_mats.part_func_mat_4_base_pairings.insert(pp_closing_loop, sum);
          let sum = sum + fe_multi;
          ss_part_func_mats.part_func_mat_4_base_pairings_accessible_on_el.insert(pp_closing_loop, sum + feature_score_sets.external_loop_accessible_basepairing_count);
          ss_part_func_mats.part_func_mat_4_base_pairings_accessible_on_mls.insert(pp_closing_loop, sum + feature_score_sets.multi_loop_basepairing_count);
        }
      }
      sum = NEG_INFINITY;
      let mut sum_2 = sum;
      for k in range_inclusive(i + T::one(), j) {
        let accessible_pp = (i, k);
        match ss_part_func_mats.part_func_mat_4_base_pairings_accessible_on_el.get(&accessible_pp) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + feature_score_sets.external_loop_accessible_baseunpairing_count * (j - k).to_f32().unwrap());
            logsumexp(&mut sum_2, ss_part_func_mats.part_func_mat_4_base_pairings_accessible_on_mls[&accessible_pp] + feature_score_sets.multi_loop_accessible_baseunpairing_count * (j - k).to_f32().unwrap());
          }, None => {},
        }
      }
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_el[long_i][long_j] = sum;
      ss_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[long_i][long_j] = sum_2;
      sum = feature_score_sets.external_loop_accessible_baseunpairing_count * sub_seq_len.to_f32().unwrap();
      for k in long_i .. long_j {
        let ss_part_func_4_rightmost_base_pairings_on_el = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_el[k][long_j];
        let part_func = ss_part_func_mats.part_func_mat[long_i][k - 1];
        logsumexp(&mut sum, part_func + ss_part_func_4_rightmost_base_pairings_on_el);
      }
      ss_part_func_mats.part_func_mat[long_i][long_j] = sum;
      sum = NEG_INFINITY;
      sum_2 = sum;
      for k in long_i .. long_j {
        let ss_part_func_4_rightmost_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[k][long_j];
        logsumexp(&mut sum, ss_part_func_4_rightmost_base_pairings_on_mls + feature_score_sets.multi_loop_accessible_baseunpairing_count * (k - long_i) as FreeEnergy);
        logsumexp(&mut sum_2, ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][k - 1] + ss_part_func_4_rightmost_base_pairings_on_mls);
      }
      ss_part_func_mats.part_func_mat_4_first_base_pairings_on_mls[long_i][long_j] = sum;
      ss_part_func_mats.part_func_mat_4_ml[long_i][long_j] = sum_2;
      logsumexp(&mut sum, sum_2);
      ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][long_j] = sum;
    }
  }
  ss_part_func_mats
}

fn get_base_pairing_prob_mat_trained<T>(seq: SeqSlice, feature_score_sets: &FeatureCountSets, trains_score_params: bool, ss_part_func_mats: &SsPartFuncMatsContra<T>, ss_free_energy_mats: &SsFreeEnergyMats<T>, allows_short_hairpins: bool, expected_feature_count_sets: &mut FeatureCountSets,) -> (SparseProbMat<T>, PartFunc)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let seq_len = seq.len();
  let ss_part_func = ss_part_func_mats.part_func_mat[1][seq_len - 2];
  let mut bpp_mat = SparseProbMat::<T>::default();
  let mut prob_mat_4_mls_1 = vec![vec![NEG_INFINITY; seq_len]; seq_len];
  let mut prob_mat_4_mls_2 = prob_mat_4_mls_1.clone();
  let short_seq_len = T::from_usize(seq_len).unwrap();
  for sub_seq_len in range_inclusive(T::from_usize(if allows_short_hairpins {2} else {MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL}).unwrap(), short_seq_len - T::from_usize(2).unwrap()).rev() {
    for i in range(T::one(), short_seq_len - sub_seq_len) {
      let j = i + sub_seq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let mut sum_1 = NEG_INFINITY;
      let mut sum_2 = sum_1;
      for k in range_inclusive(j + T::one(), short_seq_len - T::from_usize(2).unwrap()) {
        let long_k = k.to_usize().unwrap();
        let pp_closing_loop = (i, k);
        match ss_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
          Some(&part_func) => {
            let bpp = bpp_mat[&pp_closing_loop];
            let coefficient = bpp + ss_free_energy_mats.ml_closing_bp_fe_mat[&pp_closing_loop] - part_func;
            logsumexp(&mut sum_1, coefficient + ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_j + 1][long_k - 1]);
            logsumexp(&mut sum_2, coefficient + feature_score_sets.multi_loop_accessible_baseunpairing_count * (k - j - T::one()).to_f32().unwrap());
          }, None => {},
        }
      }
      prob_mat_4_mls_1[long_i][long_j] = sum_1;
      prob_mat_4_mls_2[long_i][long_j] = sum_2;
      let accessible_pp = (i, j);
      match ss_part_func_mats.part_func_mat_4_base_pairings.get(&accessible_pp) {
        Some(&part_func) => {
          let part_func_pair = (
            ss_part_func_mats.part_func_mat[1][long_i - 1],
            ss_part_func_mats.part_func_mat[long_j + 1][seq_len - 2],
          );
          let mut sum = part_func_pair.0 + part_func_pair.1 + ss_part_func_mats.part_func_mat_4_base_pairings_accessible_on_el[&accessible_pp] - ss_part_func;
          assert!(NEG_INFINITY <= sum && sum <= 0.);
          let base_pair = (seq[long_i], seq[long_j]);
          let mismatch_pair = (seq[long_j + 1], seq[long_i - 1]);
          if trains_score_params {
            // Count external loop accessible basepairings.
            logsumexp(&mut expected_feature_count_sets
              .external_loop_accessible_basepairing_count, sum);
            // Count helix ends.
            logsumexp(&mut expected_feature_count_sets.helix_end_count_mat[base_pair.1][base_pair.0], sum);
            // Count external loop terminal mismatches.
            if j < short_seq_len - T::from_usize(2).unwrap() {
              logsumexp(&mut expected_feature_count_sets.left_dangle_count_mat
                [base_pair.1][base_pair.0][seq[long_j + 1]], sum);
            }
            if i > T::one() {
              logsumexp(&mut expected_feature_count_sets.right_dangle_count_mat
                [base_pair.1][base_pair.0][seq[long_i - 1]], sum);
            }
          }
          for k in range(T::one(), i).rev() {
            let long_k = k.to_usize().unwrap();
            if long_i - long_k - 1 > CONTRA_MAX_LOOP_LEN {break;}
            for l in range_inclusive(j + T::one(), short_seq_len - T::from_usize(2).unwrap()) {
              let long_l = l.to_usize().unwrap();
              if long_l - long_j - 1 + long_i - long_k - 1 > CONTRA_MAX_LOOP_LEN {break;}
              let pp_closing_loop = (k, l);
              match ss_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
                Some(&part_func_2) => {
                  let bpp_4_2l = bpp_mat[&pp_closing_loop] + part_func - part_func_2 + ss_free_energy_mats.twoloop_fe_4d_mat[&(k, l, i, j)];
                  assert!(NEG_INFINITY <= bpp_4_2l && bpp_4_2l <= 0.);
                  logsumexp(&mut sum, bpp_4_2l);
                  if trains_score_params {
                    let loop_len_pair = (long_i - long_k - 1, long_l - long_j - 1);
                    let is_stack = loop_len_pair.0 == 0 && loop_len_pair.1 == 0;
                    let base_pair_2 = (seq[long_k], seq[long_l]);
                    let is_bulge_loop = (loop_len_pair.0 == 0 || loop_len_pair.1 == 0) && !is_stack;
                    let mismatch_pair_2 = (seq[long_k + 1], seq[long_l - 1]);
                    if is_stack {
                      // Count a stack.
                      let dict_min_stack = get_dict_min_stack(&base_pair_2, &base_pair);
                      logsumexp(&mut expected_feature_count_sets.stack_count_mat[dict_min_stack.0.0]
                        [dict_min_stack.0.1][dict_min_stack.1.0][dict_min_stack.1.1], bpp_4_2l);
                    } else {
                      if is_bulge_loop {
                        // Count a bulge loop length.
                        let bulge_loop_len = if loop_len_pair.0 == 0 {
                          loop_len_pair.1
                        } else {
                          loop_len_pair.0
                        };
                        logsumexp(&mut expected_feature_count_sets.bulge_loop_length_counts
                          [bulge_loop_len - 1], bpp_4_2l);
                        // Count a 0x1 bulge loop.
                        if bulge_loop_len == 1 {
                          let mismatch = if loop_len_pair.0 == 0 {mismatch_pair_2.1} else {mismatch_pair_2.0};
                          logsumexp(&mut expected_feature_count_sets.bulge_loop_0x1_length_counts
                            [mismatch], bpp_4_2l);
                        }
                      } else {
                        // Count an interior loop length.
                        logsumexp(&mut expected_feature_count_sets.interior_loop_length_counts
                          [loop_len_pair.0 + loop_len_pair.1 - 2], bpp_4_2l);
                        let diff = get_diff(loop_len_pair.0, loop_len_pair.1);
                        if diff == 0 {
                          logsumexp(&mut expected_feature_count_sets.interior_loop_length_counts_symm
                            [loop_len_pair.0 - 1], bpp_4_2l);
                        } else {
                          logsumexp(&mut expected_feature_count_sets.interior_loop_length_counts_asymm
                            [diff - 1], bpp_4_2l);
                        }
                        // Count a 1x1 interior loop.
                        if loop_len_pair.0 == 1 && loop_len_pair.1 == 1 {
                          let dict_min_mismatch_pair_2 = get_dict_min_mismatch_pair(&mismatch_pair_2);
                          logsumexp(&mut expected_feature_count_sets.interior_loop_1x1_length_count_mat
                            [dict_min_mismatch_pair_2.0][dict_min_mismatch_pair_2.1], bpp_4_2l);
                        }
                        // Count an explicit interior loop length pair.
                        if loop_len_pair.0 <= CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT && loop_len_pair.1 <= CONSPROB_MAX_INTERIOR_LOOP_LEN_EXPLICIT {
                          let dict_min_loop_len_pair =  get_dict_min_loop_len_pair(&loop_len_pair);
                          logsumexp(&mut expected_feature_count_sets.interior_loop_length_count_mat_explicit
                            [dict_min_loop_len_pair.0 - 1][dict_min_loop_len_pair.1 - 1], bpp_4_2l);
                        }
                      }
                      // Count helix ends.
                      logsumexp(&mut expected_feature_count_sets.helix_end_count_mat[base_pair.1][base_pair.0], bpp_4_2l);
                      logsumexp(&mut expected_feature_count_sets.helix_end_count_mat[base_pair_2.0][base_pair_2.1], bpp_4_2l);
                      // Count 2-loop terminal mismatches.
                      logsumexp(&mut expected_feature_count_sets.terminal_mismatch_count_mat
                        [base_pair_2.0][base_pair_2.1][mismatch_pair_2.0]
                        [mismatch_pair_2.1], bpp_4_2l);
                      logsumexp(&mut expected_feature_count_sets.terminal_mismatch_count_mat
                        [base_pair.1][base_pair.0][mismatch_pair.0][mismatch_pair.1],
                        bpp_4_2l);
                    }
                  }
                }, None => {},
              }
            }
          }
          let coefficient = ss_part_func_mats.part_func_mat_4_base_pairings_accessible_on_mls[&accessible_pp];
          let mut bpp_4_ml = NEG_INFINITY;
          for k in 1 .. long_i {
            let ss_part_func_4_at_least_1_base_pairings_on_mls = ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_i - 1];
            logsumexp(&mut bpp_4_ml, coefficient + prob_mat_4_mls_2[k][long_j] + ss_part_func_4_at_least_1_base_pairings_on_mls);
            let prob_4_mls = prob_mat_4_mls_1[k][long_j];
            logsumexp(&mut bpp_4_ml, coefficient + prob_4_mls + feature_score_sets.multi_loop_accessible_baseunpairing_count * (long_i - k - 1) as FreeEnergy);
            logsumexp(&mut bpp_4_ml, coefficient + prob_4_mls + ss_part_func_4_at_least_1_base_pairings_on_mls);
          }
          assert!(NEG_INFINITY <= bpp_4_ml && bpp_4_ml <= 0.);
          logsumexp(&mut sum, bpp_4_ml);
          if trains_score_params && bpp_4_ml > NEG_INFINITY {
            // Count multi-loop terminal mismatches.
            logsumexp(&mut expected_feature_count_sets.left_dangle_count_mat
              [base_pair.1][base_pair.0][mismatch_pair.0],
              bpp_4_ml);
            logsumexp(&mut expected_feature_count_sets.right_dangle_count_mat
              [base_pair.1][base_pair.0][mismatch_pair.1],
              bpp_4_ml);
            // Count helix ends.
            logsumexp(&mut expected_feature_count_sets.helix_end_count_mat
              [base_pair.1][base_pair.0],
              bpp_4_ml);
            // Count multi-loop closings.
            logsumexp(&mut expected_feature_count_sets.multi_loop_base_count, bpp_4_ml);
            // Count multi-loop closing basepairings.
            logsumexp(&mut expected_feature_count_sets.multi_loop_basepairing_count, bpp_4_ml);
            // Count multi-loop accessible basepairings.
            logsumexp(&mut expected_feature_count_sets
              .multi_loop_basepairing_count, bpp_4_ml);
          }
          // debug_assert!(NEG_INFINITY <= sum && sum <= 0.);
          assert!(NEG_INFINITY <= sum && sum <= 0.);
          bpp_mat.insert(accessible_pp, sum);
          if trains_score_params {
            // Count base pairs.
            let dict_min_base_pair = get_dict_min_base_pair(&base_pair);
            logsumexp(&mut expected_feature_count_sets.base_pair_count_mat[dict_min_base_pair.0][dict_min_base_pair.1], sum);
            let mismatch_pair = (seq[long_i + 1], seq[long_j - 1]);
            if sub_seq_len.to_usize().unwrap() - 2 <= CONSPROB_MAX_HAIRPIN_LOOP_LEN {
              let hairpin_loop_score = ss_free_energy_mats.hl_fe_mat[&(i, j)];
              let bpp_4_hl = sum - part_func + hairpin_loop_score;
              assert!(NEG_INFINITY <= bpp_4_hl && bpp_4_hl <= 0.);
              logsumexp(&mut expected_feature_count_sets.hairpin_loop_length_counts[long_j - long_i - 1], bpp_4_hl);
              logsumexp(&mut expected_feature_count_sets.terminal_mismatch_count_mat
                [base_pair.0][base_pair.1][mismatch_pair.0]
                [mismatch_pair.1], bpp_4_hl);
              // Count helix ends.
              logsumexp(&mut expected_feature_count_sets.helix_end_count_mat[base_pair.0][base_pair.1], bpp_4_hl);
            }
            let part_func_4_ml = ss_part_func_mats.part_func_mat_4_ml[long_i + 1][long_j - 1];
            let ml_closing_bp_fe = ss_free_energy_mats.ml_closing_bp_fe_mat[&accessible_pp];
            let coeff_ml = sum + ml_closing_bp_fe - part_func;
            let bpp_4_ml = coeff_ml + part_func_4_ml;
            assert!(NEG_INFINITY <= bpp_4_ml && bpp_4_ml <= 0.);
            // Count multi-loop terminal mismatches.
            logsumexp(&mut expected_feature_count_sets.left_dangle_count_mat
              [base_pair.0][base_pair.1][mismatch_pair.0],
              bpp_4_ml);
            logsumexp(&mut expected_feature_count_sets.right_dangle_count_mat
              [base_pair.0][base_pair.1][mismatch_pair.1],
              bpp_4_ml);
            // Count helix ends.
            logsumexp(&mut expected_feature_count_sets.helix_end_count_mat
              [base_pair.0][base_pair.1], bpp_4_ml);
            for k in long_i + 1 .. long_j {
              let coeff = coeff_ml + feature_score_sets.multi_loop_accessible_baseunpairing_count;
              let part_func_pair = (
                feature_score_sets.multi_loop_accessible_baseunpairing_count * (k - long_i - 1).to_f32().unwrap(),
                ss_part_func_mats.part_func_mat_4_ml[k + 1][long_j - 1],
              );
              let mut prob_4_ml = coeff + part_func_pair.0 + part_func_pair.1;
              let part_func_pair = (
                ss_part_func_mats.part_func_mat_4_first_base_pairings_on_mls[long_i + 1][k - 1],
                ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_j - 1],
              );
              let term = coeff + part_func_pair.0 + part_func_pair.1;
              logsumexp(&mut prob_4_ml, term);
              let part_func_pair = (
                ss_part_func_mats.part_func_mat_4_ml[long_i + 1][k - 1],
                ss_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_j - 1],
              );
              let term = coeff + part_func_pair.0 + part_func_pair.1;
              logsumexp(&mut prob_4_ml, term);
              let term = coeff + part_func_pair.0 + feature_score_sets.multi_loop_accessible_baseunpairing_count * (long_j - k - 1).to_f32().unwrap();
              logsumexp(&mut prob_4_ml, term);
              assert!(NEG_INFINITY <= prob_4_ml && prob_4_ml <= 0.);
              // Count multi-loop accessible unpairings.
              logsumexp(&mut expected_feature_count_sets
                .multi_loop_accessible_baseunpairing_count, prob_4_ml);
            }
          }
        }, None => {},
      }
    }
  }
  if trains_score_params {
    for i in 1 .. seq_len - 1 {
      let part_func_pair = (
        ss_part_func_mats.part_func_mat[1][i - 1],
        ss_part_func_mats.part_func_mat[i + 1][seq_len - 2],
      );
      let prob_4_el = part_func_pair.0 + feature_score_sets.external_loop_accessible_baseunpairing_count + part_func_pair.1 - ss_part_func;
      // assert!(NEG_INFINITY <= prob_4_el && prob_4_el <= 0.);
      // Count external loop accessible unpairings.
      logsumexp(&mut expected_feature_count_sets
        .external_loop_accessible_baseunpairing_count, prob_4_el);
    }
    for count in expected_feature_count_sets.hairpin_loop_length_counts.iter_mut() {
      *count = expf(*count);
    }
    for count in expected_feature_count_sets.bulge_loop_length_counts.iter_mut() {
      *count = expf(*count);
    }
    for count in expected_feature_count_sets.interior_loop_length_counts.iter_mut() {
      *count = expf(*count);
    }
    for count in expected_feature_count_sets.interior_loop_length_counts_symm.iter_mut() {
      *count = expf(*count);
    }
    for count in expected_feature_count_sets.interior_loop_length_counts_asymm.iter_mut() {
      *count = expf(*count);
    }
    for count_3d_mat in expected_feature_count_sets.stack_count_mat.iter_mut() {
      for count_2d_mat in count_3d_mat.iter_mut() {
        for counts in count_2d_mat.iter_mut() {
          for count in counts.iter_mut() {
            *count = expf(*count);
          }
        }
      }
    }
    for count_3d_mat in expected_feature_count_sets.terminal_mismatch_count_mat.iter_mut() {
      for count_2d_mat in count_3d_mat.iter_mut() {
        for counts in count_2d_mat.iter_mut() {
          for count in counts.iter_mut() {
            *count = expf(*count);
          }
        }
      }
    }
    for count_2d_mat in expected_feature_count_sets.left_dangle_count_mat.iter_mut() {
      for counts in count_2d_mat.iter_mut() {
        for count in counts.iter_mut() {
          *count = expf(*count);
        }
      }
    }
    for count_2d_mat in expected_feature_count_sets.right_dangle_count_mat.iter_mut() {
      for counts in count_2d_mat.iter_mut() {
        for count in counts.iter_mut() {
          *count = expf(*count);
        }
      }
    }
    for counts in expected_feature_count_sets.interior_loop_length_count_mat_explicit.iter_mut() {
      for count in counts.iter_mut() {
        *count = expf(*count);
      }
    }
    for count in expected_feature_count_sets.bulge_loop_0x1_length_counts.iter_mut() {
      *count = expf(*count);
    }
    for counts in expected_feature_count_sets.interior_loop_1x1_length_count_mat.iter_mut() {
      for count in counts.iter_mut() {
        *count = expf(*count);
      }
    }
    for counts in expected_feature_count_sets.helix_end_count_mat.iter_mut() {
      for count in counts.iter_mut() {
        *count = expf(*count);
      }
    }
    for counts in expected_feature_count_sets.base_pair_count_mat.iter_mut() {
      for count in counts.iter_mut() {
        *count = expf(*count);
      }
    }
    expected_feature_count_sets.multi_loop_base_count = expf(expected_feature_count_sets.multi_loop_base_count);
    expected_feature_count_sets.multi_loop_basepairing_count = expf(expected_feature_count_sets.multi_loop_basepairing_count);
    expected_feature_count_sets.multi_loop_accessible_baseunpairing_count = expf(expected_feature_count_sets.multi_loop_accessible_baseunpairing_count);
    expected_feature_count_sets.external_loop_accessible_basepairing_count = expf(expected_feature_count_sets.external_loop_accessible_basepairing_count);
    expected_feature_count_sets.external_loop_accessible_baseunpairing_count = expf(expected_feature_count_sets.external_loop_accessible_baseunpairing_count);
  }
  bpp_mat = bpp_mat.iter().map(|(pos_pair, &bpp)| {(*pos_pair, expf(bpp))}).collect();
  (bpp_mat, ss_part_func)
}

pub fn get_hl_fe_trained(feature_score_sets: &FeatureCountSets, seq: SeqSlice, pp_closing_loop: &(usize, usize)) -> FreeEnergy {
  let hl_len = pp_closing_loop.1 - pp_closing_loop.0 - 1;
  feature_score_sets.hairpin_loop_length_counts_cumulative[hl_len]
    + get_junction_fe_single_trained(feature_score_sets, seq, pp_closing_loop)
}

pub fn get_2_loop_fe_trained(
  feature_score_sets: &FeatureCountSets,
  seq: SeqSlice,
  pp_closing_loop: &(usize, usize),
  accessible_pp: &(usize, usize),
) -> FeatureCount {
    let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  let fe = if pp_closing_loop.0 + 1 == accessible_pp.0 && pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_stack_fe_trained(feature_score_sets, seq, pp_closing_loop, accessible_pp)
  } else if pp_closing_loop.0 + 1 == accessible_pp.0 || pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_bl_fe_trained(feature_score_sets, seq, pp_closing_loop, accessible_pp)
  } else {
    get_il_fe_trained(feature_score_sets, seq, pp_closing_loop, accessible_pp)
  };
  fe + feature_score_sets.base_pair_count_mat[accessible_bp.0][accessible_bp.1]
}

pub fn get_stack_fe_trained(feature_score_sets: &FeatureCountSets, seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  feature_score_sets.stack_count_mat[bp_closing_loop.0][bp_closing_loop.1][accessible_bp.0][accessible_bp.1]
}

pub fn get_bl_fe_trained(feature_score_sets: &FeatureCountSets, seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let bl_len = accessible_pp.0 - pp_closing_loop.0 + pp_closing_loop.1 - accessible_pp.1 - 2;
  let fe = if bl_len == 1 {
    feature_score_sets.bulge_loop_0x1_length_counts[if accessible_pp.0 - pp_closing_loop.0 - 1 == 1 {
      seq[pp_closing_loop.0 + 1]
    } else {
      seq[pp_closing_loop.1 - 1]
    }]
  } else {0.};
  fe + feature_score_sets.bulge_loop_length_counts_cumulative[bl_len - 1]
    + get_junction_fe_single_trained(feature_score_sets, seq, pp_closing_loop)
    + get_junction_fe_single_trained(feature_score_sets, seq, &(accessible_pp.1, accessible_pp.0))
}

pub fn get_il_fe_trained(feature_score_sets: &FeatureCountSets, seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let pair_of_nums_of_unpaired_bases = (accessible_pp.0 - pp_closing_loop.0 - 1, pp_closing_loop.1 - accessible_pp.1 - 1);
  let il_len = pair_of_nums_of_unpaired_bases.0 + pair_of_nums_of_unpaired_bases.1;
  let fe = if pair_of_nums_of_unpaired_bases.0 == pair_of_nums_of_unpaired_bases.1 {
    let fe_3 = if il_len == 2 {
      feature_score_sets.interior_loop_1x1_length_count_mat[seq[pp_closing_loop.0 + 1]][seq[pp_closing_loop.1 - 1]]
    } else {0.};
    fe_3 + feature_score_sets.interior_loop_length_counts_symm_cumulative[pair_of_nums_of_unpaired_bases.0 - 1]
  } else {
    feature_score_sets.interior_loop_length_counts_asymm_cumulative[get_abs_diff(pair_of_nums_of_unpaired_bases.0, pair_of_nums_of_unpaired_bases.1) - 1]
  };
  let fe_2 = if pair_of_nums_of_unpaired_bases.0 <= 4 && pair_of_nums_of_unpaired_bases.1 <= 4 {
    feature_score_sets.interior_loop_length_count_mat_explicit[pair_of_nums_of_unpaired_bases.0 - 1][pair_of_nums_of_unpaired_bases.1 - 1]
  } else {0.};
  fe + fe_2 + feature_score_sets.interior_loop_length_counts_cumulative[il_len - 2]
    + get_junction_fe_single_trained(feature_score_sets, seq, pp_closing_loop)
    + get_junction_fe_single_trained(feature_score_sets, seq, &(accessible_pp.1, accessible_pp.0))
}

pub fn get_junction_fe_single_trained(feature_score_sets: &FeatureCountSets, seq: SeqSlice, pp: &(usize, usize)) -> FreeEnergy {
  let bp = (seq[pp.0], seq[pp.1]);
  get_helix_closing_fe_trained(feature_score_sets, &bp) + get_terminal_mismatch_fe_trained(feature_score_sets, &bp, &(seq[pp.0 + 1], seq[pp.1 - 1]))
}

pub fn get_helix_closing_fe_trained(feature_score_sets: &FeatureCountSets, bp: &BasePair) -> FreeEnergy {
  feature_score_sets.helix_end_count_mat[bp.0][bp.1]
}

pub fn get_terminal_mismatch_fe_trained(feature_score_sets: &FeatureCountSets, bp: &BasePair, mismatch_bp: &BasePair) -> FreeEnergy {
  feature_score_sets.terminal_mismatch_count_mat[bp.0][bp.1][mismatch_bp.0][mismatch_bp.1]
}

pub fn get_junction_fe_multi_trained(feature_score_sets: &FeatureCountSets, seq: SeqSlice, pp: &(usize, usize), seq_len: usize,) -> FreeEnergy {
  let bp = (seq[pp.0], seq[pp.1]);
  let five_prime_end = 1;
  let three_prime_end = seq_len - 2;
  get_helix_closing_fe_trained(feature_score_sets, &bp) + if pp.0 < three_prime_end && pp.1 > five_prime_end {
    feature_score_sets.left_dangle_count_mat[bp.0][bp.1][seq[pp.0 + 1]] + feature_score_sets.right_dangle_count_mat[bp.0][bp.1][seq[pp.1 - 1]]
  } else if pp.0 < three_prime_end {
    feature_score_sets.left_dangle_count_mat[bp.0][bp.1][seq[pp.0 + 1]]
  } else if pp.1 > five_prime_end {
    feature_score_sets.right_dangle_count_mat[bp.0][bp.1][seq[pp.1 - 1]]
  } else {
    0.
  }
}

pub fn get_dict_min_stack(base_pair_closing_loop: &BasePair, base_pair_accessible: &BasePair) -> (BasePair, BasePair) {
  let stack = (*base_pair_closing_loop, *base_pair_accessible);
  let inverse_stack = ((base_pair_accessible.1, base_pair_accessible.0), (base_pair_closing_loop.1, base_pair_closing_loop.0));
  if stack < inverse_stack {
    stack
  } else {
    inverse_stack
  }
}

pub fn get_dict_min_loop_len_pair(loop_len_pair: &(usize, usize)) -> (usize, usize) {
  let inverse_loop_len_pair = (loop_len_pair.1, loop_len_pair.0);
  if *loop_len_pair < inverse_loop_len_pair {
    *loop_len_pair
  } else {
    inverse_loop_len_pair
  }
}

pub fn get_dict_min_nuc_pair(nuc_pair: &(usize, usize)) -> (usize, usize) {
  get_dict_min_loop_len_pair(nuc_pair)
}

pub fn get_dict_min_mismatch_pair(mismatch_pair: &BasePair) -> BasePair {
  get_dict_min_loop_align(mismatch_pair)
}

pub fn get_dict_min_base_pair(base_pair: &BasePair) -> BasePair {
  get_dict_min_loop_align(base_pair)
}

pub fn get_num_of_multiloop_baseunpairing_nucs(pos_pair_closing_loop: &(usize, usize), pos_pairs_in_loop: &Vec<(usize, usize)>) -> usize {
  let mut count = 0;
  let mut prev_pos_pair = (0, 0);
  for pos_pair_in_loop in pos_pairs_in_loop {
    let start = if prev_pos_pair == (0, 0) {
      pos_pair_closing_loop.0 + 1
    } else {
      prev_pos_pair.1 + 1
    };
    count += pos_pair_in_loop.0 - start;
    prev_pos_pair = *pos_pair_in_loop;
  }
  count += pos_pair_closing_loop.1 - (prev_pos_pair.1 + 1);
  count
}

pub fn get_num_of_externalloop_baseunpairing_nucs(pos_pairs_in_loop: &Vec<(usize, usize)>, seq: SeqSlice) -> usize {
  let mut count = 0;
  let mut prev_pos_pair = (0, 0);
  for pos_pair_in_loop in pos_pairs_in_loop {
    let start = if prev_pos_pair == (0, 0) {
      1
    } else {
      prev_pos_pair.1 + 1
    };
    count += pos_pair_in_loop.0 - start;
    prev_pos_pair = *pos_pair_in_loop;
  }
  count += seq.len() - 1 - (prev_pos_pair.1 + 1);
  count
}

pub fn get_mismatch_pair(seq: SeqSlice, pos_pair: &(usize, usize), is_closing: bool) -> (usize, usize) {
  let mut mismatch_pair = (PSEUDO_BASE, PSEUDO_BASE);
  if is_closing {
    mismatch_pair.0 = seq[pos_pair.0 + 1];
    mismatch_pair.1 = seq[pos_pair.1 - 1];
  } else {
    mismatch_pair.0 = seq[pos_pair.0 - 1];
    mismatch_pair.1 = seq[pos_pair.1 + 1];
  }
  mismatch_pair
}

pub fn get_diff(x: usize, y: usize) -> usize {
  max(x, y) - min(x, y)
}

pub fn convert_vec_2_struct(feature_counts: &FeatureCounts, uses_cumulative_feature_counts: bool) -> FeatureCountSets {
  let mut f = FeatureCountSets::new(0.);
  let mut offset = 0;
  let len = f.hairpin_loop_length_counts.len();
  for i in 0 .. len {
    let count = feature_counts[offset + i];
    if uses_cumulative_feature_counts {
      f.hairpin_loop_length_counts_cumulative[i] = count;
    } else {
      f.hairpin_loop_length_counts[i] = count;
    }
  }
  offset += len;
  let len = f.bulge_loop_length_counts.len();
  for i in 0 .. len {
    let count = feature_counts[offset + i];
    if uses_cumulative_feature_counts {
      f.bulge_loop_length_counts_cumulative[i] = count;
    } else {
      f.bulge_loop_length_counts[i] = count;
    }
  }
  offset += len;
  let len = f.interior_loop_length_counts.len();
  for i in 0 .. len {
    let count = feature_counts[offset + i];
    if uses_cumulative_feature_counts {
      f.interior_loop_length_counts_cumulative[i] = count;
    } else {
      f.interior_loop_length_counts[i] = count;
    }
  }
  offset += len;
  let len = f.interior_loop_length_counts_symm.len();
  for i in 0 .. len {
    let count = feature_counts[offset + i];
    if uses_cumulative_feature_counts {
      f.interior_loop_length_counts_symm_cumulative[i] = count;
    } else {
      f.interior_loop_length_counts_symm[i] = count;
    }
  }
  offset += len;
  let len = f.interior_loop_length_counts_asymm.len();
  for i in 0 .. len {
    let count = feature_counts[offset + i];
    if uses_cumulative_feature_counts {
      f.interior_loop_length_counts_asymm_cumulative[i] = count;
    } else {
      f.interior_loop_length_counts_asymm[i] = count;
    }
  }
  offset += len;
  let len = f.stack_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        for l in 0 .. len {
          f.stack_count_mat[i][j][k][l] = feature_counts[offset + i * len.pow(3) + j * len.pow(2) + k * len + l];
        }
      }
    }
  }
  offset += len.pow(4);
  let len = f.terminal_mismatch_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        for l in 0 .. len {
          f.terminal_mismatch_count_mat[i][j][k][l] = feature_counts[offset + i * len.pow(3) + j * len.pow(2) + k * len + l];
        }
      }
    }
  }
  offset += len.pow(4);
  let len = f.left_dangle_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        f.left_dangle_count_mat[i][j][k] = feature_counts[offset + i * len.pow(2) + j * len + k];
      }
    }
  }
  offset += len.pow(3);
  let len = f.right_dangle_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        f.right_dangle_count_mat[i][j][k] = feature_counts[offset + i * len.pow(2) + j * len + k];
      }
    }
  }
  offset += len.pow(3);
  let len = f.helix_end_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      f.helix_end_count_mat[i][j] = feature_counts[offset + i * len + j];
    }
  }
  offset += len.pow(2);
  let len = f.base_pair_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      f.base_pair_count_mat[i][j] = feature_counts[offset + i * len + j];
    }
  }
  offset += len.pow(2);
  let len = f.interior_loop_length_count_mat_explicit.len();
  for i in 0 .. len {
    for j in 0 .. len {
      f.interior_loop_length_count_mat_explicit[i][j] = feature_counts[offset + i * len + j];
    }
  }
  offset += len.pow(2);
  let len = f.bulge_loop_0x1_length_counts.len();
  for i in 0 .. len {
    f.bulge_loop_0x1_length_counts[i] = feature_counts[offset + i];
  }
  offset += len;
  let len = f.interior_loop_1x1_length_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      f.interior_loop_1x1_length_count_mat[i][j] = feature_counts[offset + i * len + j];
    }
  }
  offset += len.pow(2);
  f.multi_loop_base_count = feature_counts[offset];
  offset += 1;
  f.multi_loop_basepairing_count = feature_counts[offset];
  offset += 1;
  f.multi_loop_accessible_baseunpairing_count = feature_counts[offset];
  offset += 1;
  f.external_loop_accessible_basepairing_count = feature_counts[offset];
  offset += 1;
  f.external_loop_accessible_baseunpairing_count = feature_counts[offset];
  offset += 1;
  assert!(offset == f.get_len());
  f
}

pub fn convert_struct_2_vec(feature_count_sets: &FeatureCountSets, uses_cumulative_feature_counts: bool) -> FeatureCounts {
  let f = feature_count_sets;
  let mut feature_counts = vec![0.; f.get_len()];
  let mut offset = 0;
  let len = f.hairpin_loop_length_counts.len();
  for i in 0 .. len {
    feature_counts[offset + i] = if uses_cumulative_feature_counts {
      f.hairpin_loop_length_counts_cumulative[i]
    } else {
      f.hairpin_loop_length_counts[i]
    };
  }
  offset += len;
  let len = f.bulge_loop_length_counts.len();
  for i in 0 .. len {
    feature_counts[offset + i] = if uses_cumulative_feature_counts {
      f.bulge_loop_length_counts_cumulative[i]
    } else {
      f.bulge_loop_length_counts[i]
    };
  }
  offset += len;
  let len = f.interior_loop_length_counts.len();
  for i in 0 .. len {
    feature_counts[offset + i] = if uses_cumulative_feature_counts {
      f.interior_loop_length_counts_cumulative[i]
    } else {
      f.interior_loop_length_counts[i]
    };
  }
  offset += len;
  let len = f.interior_loop_length_counts_symm.len();
  for i in 0 .. len {
    feature_counts[offset + i] = if uses_cumulative_feature_counts {
      f.interior_loop_length_counts_symm_cumulative[i]
    } else {
      f.interior_loop_length_counts_symm[i]
    };
  }
  offset += len;
  let len = f.interior_loop_length_counts_asymm.len();
  for i in 0 .. len {
    feature_counts[offset + i] = if uses_cumulative_feature_counts {
      f.interior_loop_length_counts_asymm_cumulative[i]
    } else {
      f.interior_loop_length_counts_asymm[i]
    };
  }
  offset += len;
  let len = f.stack_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        for l in 0 .. len {
          feature_counts[offset + i * len.pow(3) + j * len.pow(2) + k * len + l] = f.stack_count_mat[i][j][k][l];
        }
      }
    }
  }
  offset += len.pow(4);
  let len = f.terminal_mismatch_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        for l in 0 .. len {
          feature_counts[offset + i * len.pow(3) + j * len.pow(2) + k * len + l] = f.terminal_mismatch_count_mat[i][j][k][l];
        }
      }
    }
  }
  offset += len.pow(4);
  let len = f.left_dangle_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        feature_counts[offset + i * len.pow(2) + j * len + k] = f.left_dangle_count_mat[i][j][k];
      }
    }
  }
  offset += len.pow(3);
  let len = f.right_dangle_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      for k in 0 .. len {
        feature_counts[offset + i * len.pow(2) + j * len + k] = f.right_dangle_count_mat[i][j][k];
      }
    }
  }
  offset += len.pow(3);
  let len = f.helix_end_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      feature_counts[offset + i * len + j] = f.helix_end_count_mat[i][j];
    }
  }
  offset += len.pow(2);
  let len = f.base_pair_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      feature_counts[offset + i * len + j] = f.base_pair_count_mat[i][j];
    }
  }
  offset += len.pow(2);
  let len = f.interior_loop_length_count_mat_explicit.len();
  for i in 0 .. len {
    for j in 0 .. len {
      feature_counts[offset + i * len + j] = f.interior_loop_length_count_mat_explicit[i][j];
    }
  }
  offset += len.pow(2);
  let len = f.bulge_loop_0x1_length_counts.len();
  for i in 0 .. len {
    feature_counts[offset + i] = f.bulge_loop_0x1_length_counts[i];
  }
  offset += len;
  let len = f.interior_loop_1x1_length_count_mat.len();
  for i in 0 .. len {
    for j in 0 .. len {
      feature_counts[offset + i * len + j] = f.interior_loop_1x1_length_count_mat[i][j];
    }
  }
  offset += len.pow(2);
  feature_counts[offset] = f.multi_loop_base_count;
  offset += 1;
  feature_counts[offset] = f.multi_loop_basepairing_count;
  offset += 1;
  feature_counts[offset] = f.multi_loop_accessible_baseunpairing_count;
  offset += 1;
  feature_counts[offset] = f.external_loop_accessible_basepairing_count;
  offset += 1;
  feature_counts[offset] = f.external_loop_accessible_baseunpairing_count;
  offset += 1;
  assert!(offset == f.get_len());
  Array::from(feature_counts)
}

pub fn write_feature_score_sets_trained(feature_score_sets: &FeatureCountSets) {
  let mut writer_2_trained_feature_score_sets_file = BufWriter::new(File::create(TRAINED_FEATURE_SCORE_SETS_FILE_PATH).unwrap());
  let mut buf_4_writer_2_trained_feature_score_sets_file = String::from("use FeatureCountSets;\nimpl FeatureCountSets {\npub fn load_trained_score_params() -> FeatureCountSets {\nFeatureCountSets {\nhairpin_loop_length_counts: ");
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nbulge_loop_length_counts: ", &feature_score_sets.hairpin_loop_length_counts));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_length_counts: ", &feature_score_sets.bulge_loop_length_counts));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_length_counts_symm: ", &feature_score_sets.interior_loop_length_counts));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_length_counts_asymm: ", &feature_score_sets.interior_loop_length_counts_symm));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nstack_count_mat: ", &feature_score_sets.interior_loop_length_counts_asymm));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nterminal_mismatch_count_mat: ", &feature_score_sets.stack_count_mat));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nleft_dangle_count_mat: ", &feature_score_sets.terminal_mismatch_count_mat));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nright_dangle_count_mat: ", &feature_score_sets.left_dangle_count_mat));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nhelix_end_count_mat: ", &feature_score_sets.right_dangle_count_mat));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nbase_pair_count_mat: ", &feature_score_sets.helix_end_count_mat));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_length_count_mat_explicit: ", &feature_score_sets.base_pair_count_mat));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nbulge_loop_0x1_length_counts: ", &feature_score_sets.interior_loop_length_count_mat_explicit));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_1x1_length_count_mat: ", &feature_score_sets.bulge_loop_0x1_length_counts));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nmulti_loop_base_count: ", &feature_score_sets.interior_loop_1x1_length_count_mat));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nmulti_loop_basepairing_count: ", feature_score_sets.multi_loop_base_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nmulti_loop_accessible_baseunpairing_count: ", feature_score_sets.multi_loop_basepairing_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nexternal_loop_accessible_basepairing_count: ", feature_score_sets.multi_loop_accessible_baseunpairing_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nexternal_loop_accessible_baseunpairing_count: ", feature_score_sets.external_loop_accessible_basepairing_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nhairpin_loop_length_counts_cumulative: ", feature_score_sets.external_loop_accessible_baseunpairing_count));
  /* buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nmatch_2_match_count: ", feature_score_sets.external_loop_accessible_baseunpairing_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nmatch_2_insert_count: ", feature_score_sets.match_2_match_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninsert_extend_count: ", feature_score_sets.match_2_insert_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninsert_switch_count: ", feature_score_sets.insert_extend_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninit_match_count: ", feature_score_sets.insert_switch_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninit_insert_count: ", feature_score_sets.init_match_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninsert_counts: ", feature_score_sets.init_insert_count));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nalign_count_mat: ", feature_score_sets.insert_counts));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nhairpin_loop_length_counts_cumulative: ", feature_score_sets.align_count_mat)); */
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\nbulge_loop_length_counts_cumulative: ", &feature_score_sets.hairpin_loop_length_counts_cumulative));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_length_counts_cumulative: ", &feature_score_sets.bulge_loop_length_counts_cumulative));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_length_counts_symm_cumulative: ", &feature_score_sets.interior_loop_length_counts_cumulative));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},\ninterior_loop_length_counts_asymm_cumulative: ", &feature_score_sets.interior_loop_length_counts_symm_cumulative));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&format!("{:?},", &feature_score_sets.interior_loop_length_counts_asymm_cumulative));
  buf_4_writer_2_trained_feature_score_sets_file.push_str(&String::from("\n}\n}\n}"));
  let _ = writer_2_trained_feature_score_sets_file.write_all(buf_4_writer_2_trained_feature_score_sets_file.as_bytes());
}

pub fn write_costs(costs: &Probs, output_file_path: &Path) {
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  let mut buf_4_writer_2_output_file = String::new();
  for cost in costs {
    buf_4_writer_2_output_file.push_str(&format!("{}\n", cost));
  }
  let _ = writer_2_output_file.write_all(buf_4_writer_2_output_file.as_bytes());
}
