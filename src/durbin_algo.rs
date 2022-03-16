use utils::*;
use compiled_seq_align_params::*;

pub struct SaPartFuncMats {
  pub forward_part_func_mat_match: PartFuncMat,
  pub forward_part_func_mat_insert: PartFuncMat,
  pub forward_part_func_mat_del: PartFuncMat,
  pub backward_part_func_mat_match: PartFuncMat,
  pub backward_part_func_mat_insert: PartFuncMat,
  pub backward_part_func_mat_del: PartFuncMat,
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

pub fn durbin_algo(seq_pair: &SeqPair) -> ProbMat {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let sa_part_func_mats = get_sa_part_func_mats(seq_pair);
  get_match_prob_mat(&sa_part_func_mats, &seq_len_pair)
}

pub fn get_sa_part_func_mats(seq_pair: &SeqPair) -> SaPartFuncMats {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let mut sa_part_func_mats = SaPartFuncMats::new(&seq_len_pair);
  for i in 0 .. seq_len_pair.0 - 1 {
    for j in 0 .. seq_len_pair.1 - 1 {
      if i == 0 && j == 0 {
        sa_part_func_mats.forward_part_func_mat_match[i][j] = 0.;
        continue;
      }
      if i > 0 && j > 0 {
        let mut sum = NEG_INFINITY;
        let char_pair = (seq_pair.0[i], seq_pair.1[j]);
        let match_score = MATCH_SCORE_MAT[char_pair.0][char_pair.1];
        let is_begin = (i - 1, j - 1) == (0, 0);
        let term = sa_part_func_mats.forward_part_func_mat_match[i - 1][j - 1] + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_insert[i - 1][j - 1] + MATCH_2_INSERT_SCORE;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_del[i - 1][j - 1] + MATCH_2_INSERT_SCORE;
        logsumexp(&mut sum, term);
        sa_part_func_mats.forward_part_func_mat_match[i][j] = sum + match_score;
      }
      if i > 0 {
        let chr = seq_pair.0[i];
        let insert_score = INSERT_SCORES[chr];
        let is_begin = (i - 1, j) == (0, 0);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.forward_part_func_mat_match[i - 1][j] + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_insert[i - 1][j] + INSERT_EXTEND_SCORE;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_del[i - 1][j] + INSERT_SWITCH_SCORE;
        logsumexp(&mut sum, term);
        sa_part_func_mats.forward_part_func_mat_insert[i][j] = sum + insert_score;
      }
      if j > 0 {
        let chr = seq_pair.1[j];
        let insert_score = INSERT_SCORES[chr];
        let is_begin = (i, j - 1) == (0, 0);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.forward_part_func_mat_match[i][j - 1] + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_insert[i][j - 1] + INSERT_SWITCH_SCORE;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.forward_part_func_mat_del[i][j - 1] + INSERT_EXTEND_SCORE;
        logsumexp(&mut sum, term);
        sa_part_func_mats.forward_part_func_mat_del[i][j] = sum + insert_score;
      }
    }
  }
  for i in (1 .. seq_len_pair.0).rev() {
    for j in (1 .. seq_len_pair.1).rev() {
      if i == seq_len_pair.0 - 1 && j == seq_len_pair.1 - 1 {
        sa_part_func_mats.backward_part_func_mat_match[i][j] = 0.;
        continue;
      }
      if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 {
        let mut sum = NEG_INFINITY;
        let char_pair = (seq_pair.0[i], seq_pair.1[j]);
        let match_score = MATCH_SCORE_MAT[char_pair.0][char_pair.1];
        let is_end = (i + 1, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let term = sa_part_func_mats.backward_part_func_mat_match[i + 1][j + 1] + if is_end {0.} else {MATCH_2_MATCH_SCORE};
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_insert[i + 1][j + 1] + MATCH_2_INSERT_SCORE;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_del[i + 1][j + 1] + MATCH_2_INSERT_SCORE;
        logsumexp(&mut sum, term);
        sa_part_func_mats.backward_part_func_mat_match[i][j] = sum + match_score;
      }
      if i < seq_len_pair.0 - 1 {
        let chr = seq_pair.0[i];
        let insert_score = INSERT_SCORES[chr];
        let is_end = (i + 1, j) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.backward_part_func_mat_match[i + 1][j] + if is_end {0.} else {MATCH_2_INSERT_SCORE};
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_insert[i + 1][j] + INSERT_EXTEND_SCORE;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_del[i + 1][j] + INSERT_SWITCH_SCORE;
        logsumexp(&mut sum, term);
        sa_part_func_mats.backward_part_func_mat_insert[i][j] = sum + insert_score;
      }
      if j < seq_len_pair.1 - 1 {
        let chr = seq_pair.1[j];
        let insert_score = INSERT_SCORES[chr];
        let is_end = (i, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let mut sum = NEG_INFINITY;
        let term = sa_part_func_mats.backward_part_func_mat_match[i][j + 1] + if is_end {0.} else {MATCH_2_INSERT_SCORE};
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_insert[i][j + 1] + INSERT_SWITCH_SCORE;
        logsumexp(&mut sum, term);
        let term = sa_part_func_mats.backward_part_func_mat_del[i][j + 1] + INSERT_EXTEND_SCORE;
        logsumexp(&mut sum, term);
        sa_part_func_mats.backward_part_func_mat_del[i][j] = sum + insert_score;
      }
    }
  }
  sa_part_func_mats
}

fn get_match_prob_mat(sa_part_func_mats: &SaPartFuncMats, seq_len_pair: &(usize, usize)) -> ProbMat {
  let mut match_prob_mat = vec![vec![0.; seq_len_pair.1]; seq_len_pair.0];
  let mut part_func = sa_part_func_mats.forward_part_func_mat_match[seq_len_pair.0 - 2][seq_len_pair.1 - 2];
  logsumexp(&mut part_func, sa_part_func_mats.forward_part_func_mat_insert[seq_len_pair.0 - 2][seq_len_pair.1 - 2]);
  logsumexp(&mut part_func, sa_part_func_mats.forward_part_func_mat_del[seq_len_pair.0 - 2][seq_len_pair.1 - 2]);
  for i in 1 .. seq_len_pair.0 - 1 {
    for j in 1 .. seq_len_pair.1 - 1 {
      let mut sum = NEG_INFINITY;
      let forward_part_func = sa_part_func_mats.forward_part_func_mat_match[i][j];
      let is_end = (i + 1, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
      let term = if is_end {0.} else {MATCH_2_MATCH_SCORE} + sa_part_func_mats.backward_part_func_mat_match[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let term = MATCH_2_INSERT_SCORE + sa_part_func_mats.backward_part_func_mat_insert[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let term = MATCH_2_INSERT_SCORE + sa_part_func_mats.backward_part_func_mat_del[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let match_prob = expf(forward_part_func + sum - part_func);
      debug_assert!(0. <= match_prob && match_prob <= 1.);
      match_prob_mat[i][j] = match_prob;
    }
  }
  match_prob_mat
}
