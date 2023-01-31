use compiled_seq_align_params::*;
use utils::*;

#[derive(Clone, Debug)]
pub struct AlignFeatureCountSets {
  pub match_2_match_count: FeatureCount,
  pub match_2_insert_count: FeatureCount,
  pub insert_extend_count: FeatureCount,
  pub insert_switch_count: FeatureCount,
  pub init_match_count: FeatureCount,
  pub init_insert_count: FeatureCount,
  pub insert_counts: InsertCounts,
  pub align_count_mat: AlignCountMat,
}

pub type AlignCountMat = [[FeatureCount; NUM_OF_BASES]; NUM_OF_BASES];
pub type InsertCounts = [FeatureCount; NUM_OF_BASES];

pub struct SaPartFuncMats {
  pub forward_part_func_mat_match: PartFuncMat,
  pub forward_part_func_mat_insert: PartFuncMat,
  pub forward_part_func_mat_del: PartFuncMat,
  pub backward_part_func_mat_match: PartFuncMat,
  pub backward_part_func_mat_insert: PartFuncMat,
  pub backward_part_func_mat_del: PartFuncMat,
}

impl AlignFeatureCountSets {
  pub fn new(init_val: FeatureCount) -> AlignFeatureCountSets {
    let init_vals = [init_val; NUM_OF_BASES];
    let twod_mat = [[init_val; NUM_OF_BASES]; NUM_OF_BASES];
    AlignFeatureCountSets {
      match_2_match_count: init_val,
      match_2_insert_count: init_val,
      init_match_count: init_val,
      insert_extend_count: init_val,
      insert_switch_count: init_val,
      init_insert_count: init_val,
      insert_counts: init_vals,
      align_count_mat: twod_mat,
    }
  }

  pub fn transfer(&mut self) {
    self.match_2_match_count = MATCH_2_MATCH_SCORE;
    self.match_2_insert_count = MATCH_2_INSERT_SCORE;
    self.insert_extend_count = INSERT_EXTEND_SCORE;
    self.insert_switch_count = INSERT_SWITCH_SCORE;
    self.init_match_count = INIT_MATCH_SCORE;
    self.init_insert_count = INIT_INSERT_SCORE;
    for (i, &v) in INSERT_SCORES.iter().enumerate() {
      self.insert_counts[i] = v;
    }
    for (i, vs) in MATCH_SCORE_MAT.iter().enumerate() {
      for (j, &v) in vs.iter().enumerate() {
        self.align_count_mat[i][j] = v;
      }
    }
  }
}

impl SaPartFuncMats {
  pub fn new(seq_len_pair: &(usize, usize)) -> SaPartFuncMats {
    let neg_inf_mat = vec![vec![NEG_INFINITY; seq_len_pair.1]; seq_len_pair.0];
    SaPartFuncMats {
      forward_part_func_mat_match: neg_inf_mat.clone(),
      forward_part_func_mat_insert: neg_inf_mat.clone(),
      forward_part_func_mat_del: neg_inf_mat.clone(),
      backward_part_func_mat_match: neg_inf_mat.clone(),
      backward_part_func_mat_insert: neg_inf_mat.clone(),
      backward_part_func_mat_del: neg_inf_mat,
    }
  }
}

pub fn durbin_algo(
  seq_pair: &SeqPair,
  align_feature_score_sets: &AlignFeatureCountSets,
) -> ProbMat {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let sa_part_func_mats = get_sa_part_func_mats(seq_pair, align_feature_score_sets);
  get_match_prob_mat(&sa_part_func_mats, &seq_len_pair, align_feature_score_sets)
}

pub fn get_sa_part_func_mats(
  seq_pair: &SeqPair,
  align_feature_score_sets: &AlignFeatureCountSets,
) -> SaPartFuncMats {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let mut sa_part_func_mats = SaPartFuncMats::new(&seq_len_pair);
  for i in 0..seq_len_pair.0 - 1 {
    for j in 0..seq_len_pair.1 - 1 {
      if i == 0 && j == 0 {
        sa_part_func_mats.forward_part_func_mat_match[i][j] = 0.;
        continue;
      }
      if i > 0 && j > 0 {
        let mut sum = NEG_INFINITY;
        let char_pair = (seq_pair.0[i], seq_pair.1[j]);
        let match_score = align_feature_score_sets.align_count_mat[char_pair.0][char_pair.1];
        let is_begin = (i - 1, j - 1) == (0, 0);
        let term = sa_part_func_mats.forward_part_func_mat_match[i - 1][j - 1]
          + if is_begin {
            align_feature_score_sets.init_match_count
          } else {
            align_feature_score_sets.match_2_match_count
          };
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_insert[i - 1][j - 1]
          + align_feature_score_sets.match_2_insert_count;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_del[i - 1][j - 1]
          + align_feature_score_sets.match_2_insert_count;
        logsumexp(&mut sum, term);
        sa_part_func_mats.forward_part_func_mat_match[i][j] = sum + match_score;
      }
      if i > 0 {
        let chr = seq_pair.0[i];
        let insert_score = align_feature_score_sets.insert_counts[chr];
        let is_begin = (i - 1, j) == (0, 0);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.forward_part_func_mat_match[i - 1][j]
          + if is_begin {
            align_feature_score_sets.init_insert_count
          } else {
            align_feature_score_sets.match_2_insert_count
          };
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_insert[i - 1][j]
          + align_feature_score_sets.insert_extend_count;
        logsumexp(&mut sum, term);
        sa_part_func_mats.forward_part_func_mat_insert[i][j] = sum + insert_score;
      }
      if j > 0 {
        let chr = seq_pair.1[j];
        let insert_score = align_feature_score_sets.insert_counts[chr];
        let is_begin = (i, j - 1) == (0, 0);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.forward_part_func_mat_match[i][j - 1]
          + if is_begin {
            align_feature_score_sets.init_insert_count
          } else {
            align_feature_score_sets.match_2_insert_count
          };
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_del[i][j - 1]
          + align_feature_score_sets.insert_extend_count;
        logsumexp(&mut sum, term);
        sa_part_func_mats.forward_part_func_mat_del[i][j] = sum + insert_score;
      }
    }
  }
  for i in (1..seq_len_pair.0).rev() {
    for j in (1..seq_len_pair.1).rev() {
      if i == seq_len_pair.0 - 1 && j == seq_len_pair.1 - 1 {
        sa_part_func_mats.backward_part_func_mat_match[i][j] = 0.;
        continue;
      }
      if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 {
        let mut sum = NEG_INFINITY;
        let char_pair = (seq_pair.0[i], seq_pair.1[j]);
        let match_score = align_feature_score_sets.align_count_mat[char_pair.0][char_pair.1];
        let is_end = (i + 1, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let term = sa_part_func_mats.backward_part_func_mat_match[i + 1][j + 1]
          + if is_end {
            0.
          } else {
            align_feature_score_sets.match_2_match_count
          };
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_insert[i + 1][j + 1]
          + align_feature_score_sets.match_2_insert_count;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_del[i + 1][j + 1]
          + align_feature_score_sets.match_2_insert_count;
        logsumexp(&mut sum, term);
        sa_part_func_mats.backward_part_func_mat_match[i][j] = sum + match_score;
      }
      if i < seq_len_pair.0 - 1 {
        let chr = seq_pair.0[i];
        let insert_score = align_feature_score_sets.insert_counts[chr];
        let is_end = (i + 1, j) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.backward_part_func_mat_match[i + 1][j]
          + if is_end {
            0.
          } else {
            align_feature_score_sets.match_2_insert_count
          };
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_insert[i + 1][j]
          + align_feature_score_sets.insert_extend_count;
        logsumexp(&mut sum, term);
        sa_part_func_mats.backward_part_func_mat_insert[i][j] = sum + insert_score;
      }
      if j < seq_len_pair.1 - 1 {
        let chr = seq_pair.1[j];
        let insert_score = align_feature_score_sets.insert_counts[chr];
        let is_end = (i, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.backward_part_func_mat_match[i][j + 1]
          + if is_end {
            0.
          } else {
            align_feature_score_sets.match_2_insert_count
          };
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_del[i][j + 1]
          + align_feature_score_sets.insert_extend_count;
        logsumexp(&mut sum, term);
        sa_part_func_mats.backward_part_func_mat_del[i][j] = sum + insert_score;
      }
    }
  }
  sa_part_func_mats
}

fn get_match_prob_mat(
  sa_part_func_mats: &SaPartFuncMats,
  seq_len_pair: &(usize, usize),
  align_feature_score_sets: &AlignFeatureCountSets,
) -> ProbMat {
  let mut match_prob_mat = vec![vec![0.; seq_len_pair.1]; seq_len_pair.0];
  let mut part_func =
    sa_part_func_mats.forward_part_func_mat_match[seq_len_pair.0 - 2][seq_len_pair.1 - 2];
  logsumexp(
    &mut part_func,
    sa_part_func_mats.forward_part_func_mat_insert[seq_len_pair.0 - 2][seq_len_pair.1 - 2],
  );
  logsumexp(
    &mut part_func,
    sa_part_func_mats.forward_part_func_mat_del[seq_len_pair.0 - 2][seq_len_pair.1 - 2],
  );
  for (i, vs) in match_prob_mat.iter_mut().enumerate() {
    if i == 0 || i == seq_len_pair.0 - 1 {
      continue;
    }
    for (j, v) in vs.iter_mut().enumerate() {
      if j == 0 || j == seq_len_pair.1 - 1 {
        continue;
      }
      let mut sum = NEG_INFINITY;
      let forward_part_func = sa_part_func_mats.forward_part_func_mat_match[i][j];
      let is_end = (i + 1, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
      let term = if is_end {
        0.
      } else {
        align_feature_score_sets.match_2_match_count
      } + sa_part_func_mats.backward_part_func_mat_match[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let term = align_feature_score_sets.match_2_insert_count
        + sa_part_func_mats.backward_part_func_mat_insert[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let term = align_feature_score_sets.match_2_insert_count
        + sa_part_func_mats.backward_part_func_mat_del[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let match_prob = expf(forward_part_func + sum - part_func);
      *v = match_prob;
    }
  }
  match_prob_mat
}
