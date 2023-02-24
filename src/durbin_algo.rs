use compiled_align_scores::*;
use utils::*;

#[derive(Clone, Debug)]
pub struct AlignScores {
  pub match2match_score: Score,
  pub match2insert_score: Score,
  pub insert_extend_score: Score,
  pub insert_switch_score: Score,
  pub init_match_score: Score,
  pub init_insert_score: Score,
  pub insert_scores: InsertScores,
  pub match_scores: MatchScores,
}

pub struct AlignSums {
  pub forward_sums_match: SumMat,
  pub forward_sums_insert: SumMat,
  pub forward_sums_del: SumMat,
  pub backward_sums_match: SumMat,
  pub backward_sums_insert: SumMat,
  pub backward_sums_del: SumMat,
}

impl AlignScores {
  pub fn new(init_val: Score) -> AlignScores {
    let init_vals = [init_val; NUM_BASES];
    let mat = [[init_val; NUM_BASES]; NUM_BASES];
    AlignScores {
      match2match_score: init_val,
      match2insert_score: init_val,
      init_match_score: init_val,
      insert_extend_score: init_val,
      insert_switch_score: init_val,
      init_insert_score: init_val,
      insert_scores: init_vals,
      match_scores: mat,
    }
  }

  pub fn transfer(&mut self) {
    self.match2match_score = MATCH2MATCH_SCORE;
    self.match2insert_score = MATCH2INSERT_SCORE;
    self.insert_extend_score = INSERT_EXTEND_SCORE;
    self.insert_switch_score = INSERT_SWITCH_SCORE;
    self.init_match_score = INIT_MATCH_SCORE;
    self.init_insert_score = INIT_INSERT_SCORE;
    for (i, &v) in INSERT_SCORES.iter().enumerate() {
      self.insert_scores[i] = v;
    }
    for (i, x) in MATCH_SCORES.iter().enumerate() {
      for (j, &x) in x.iter().enumerate() {
        self.match_scores[i][j] = x;
      }
    }
  }
}

impl AlignSums {
  pub fn new(seq_len_pair: &(usize, usize)) -> AlignSums {
    let neg_infs = vec![vec![NEG_INFINITY; seq_len_pair.1]; seq_len_pair.0];
    AlignSums {
      forward_sums_match: neg_infs.clone(),
      forward_sums_insert: neg_infs.clone(),
      forward_sums_del: neg_infs.clone(),
      backward_sums_match: neg_infs.clone(),
      backward_sums_insert: neg_infs.clone(),
      backward_sums_del: neg_infs,
    }
  }
}

pub fn durbin_algo(seq_pair: &SeqPair, align_scores: &AlignScores) -> ProbMat {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let align_sums = get_align_sums(seq_pair, align_scores);
  get_match_probs(&align_sums, &seq_len_pair, align_scores)
}

pub fn get_align_sums(seq_pair: &SeqPair, align_scores: &AlignScores) -> AlignSums {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let mut align_sums = AlignSums::new(&seq_len_pair);
  for i in 0..seq_len_pair.0 - 1 {
    for j in 0..seq_len_pair.1 - 1 {
      if i == 0 && j == 0 {
        align_sums.forward_sums_match[i][j] = 0.;
        continue;
      }
      if i > 0 && j > 0 {
        let mut sum = NEG_INFINITY;
        let basepair = (seq_pair.0[i], seq_pair.1[j]);
        let match_score = align_scores.match_scores[basepair.0][basepair.1];
        let begins_sum = (i - 1, j - 1) == (0, 0);
        let term = align_sums.forward_sums_match[i - 1][j - 1]
          + if begins_sum {
            align_scores.init_match_score
          } else {
            align_scores.match2match_score
          };
        logsumexp(&mut sum, term);
        let term = align_sums.forward_sums_insert[i - 1][j - 1] + align_scores.match2insert_score;
        logsumexp(&mut sum, term);
        let term = align_sums.forward_sums_del[i - 1][j - 1] + align_scores.match2insert_score;
        logsumexp(&mut sum, term);
        align_sums.forward_sums_match[i][j] = sum + match_score;
      }
      if i > 0 {
        let base = seq_pair.0[i];
        let insert_score = align_scores.insert_scores[base];
        let begins_sum = (i - 1, j) == (0, 0);
        let mut sum = NEG_INFINITY;
        let term = align_sums.forward_sums_match[i - 1][j]
          + if begins_sum {
            align_scores.init_insert_score
          } else {
            align_scores.match2insert_score
          };
        logsumexp(&mut sum, term);
        let term = align_sums.forward_sums_insert[i - 1][j] + align_scores.insert_extend_score;
        logsumexp(&mut sum, term);
        align_sums.forward_sums_insert[i][j] = sum + insert_score;
      }
      if j > 0 {
        let base = seq_pair.1[j];
        let insert_score = align_scores.insert_scores[base];
        let begins_sum = (i, j - 1) == (0, 0);
        let mut sum = NEG_INFINITY;
        let term = align_sums.forward_sums_match[i][j - 1]
          + if begins_sum {
            align_scores.init_insert_score
          } else {
            align_scores.match2insert_score
          };
        logsumexp(&mut sum, term);
        let term = align_sums.forward_sums_del[i][j - 1] + align_scores.insert_extend_score;
        logsumexp(&mut sum, term);
        align_sums.forward_sums_del[i][j] = sum + insert_score;
      }
    }
  }
  for i in (1..seq_len_pair.0).rev() {
    for j in (1..seq_len_pair.1).rev() {
      if i == seq_len_pair.0 - 1 && j == seq_len_pair.1 - 1 {
        align_sums.backward_sums_match[i][j] = 0.;
        continue;
      }
      if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 {
        let mut sum = NEG_INFINITY;
        let base_pair = (seq_pair.0[i], seq_pair.1[j]);
        let match_score = align_scores.match_scores[base_pair.0][base_pair.1];
        let ends_sum = (i + 1, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let term = align_sums.backward_sums_match[i + 1][j + 1]
          + if ends_sum {
            0.
          } else {
            align_scores.match2match_score
          };
        logsumexp(&mut sum, term);
        let term = align_sums.backward_sums_insert[i + 1][j + 1] + align_scores.match2insert_score;
        logsumexp(&mut sum, term);
        let term = align_sums.backward_sums_del[i + 1][j + 1] + align_scores.match2insert_score;
        logsumexp(&mut sum, term);
        align_sums.backward_sums_match[i][j] = sum + match_score;
      }
      if i < seq_len_pair.0 - 1 {
        let base = seq_pair.0[i];
        let insert_score = align_scores.insert_scores[base];
        let ends_sum = (i + 1, j) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let mut sum = NEG_INFINITY;
        let term = align_sums.backward_sums_match[i + 1][j]
          + if ends_sum {
            0.
          } else {
            align_scores.match2insert_score
          };
        logsumexp(&mut sum, term);
        let term = align_sums.backward_sums_insert[i + 1][j] + align_scores.insert_extend_score;
        logsumexp(&mut sum, term);
        align_sums.backward_sums_insert[i][j] = sum + insert_score;
      }
      if j < seq_len_pair.1 - 1 {
        let base = seq_pair.1[j];
        let insert_score = align_scores.insert_scores[base];
        let ends_sum = (i, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
        let mut sum = NEG_INFINITY;
        let term = align_sums.backward_sums_match[i][j + 1]
          + if ends_sum {
            0.
          } else {
            align_scores.match2insert_score
          };
        logsumexp(&mut sum, term);
        let term = align_sums.backward_sums_del[i][j + 1] + align_scores.insert_extend_score;
        logsumexp(&mut sum, term);
        align_sums.backward_sums_del[i][j] = sum + insert_score;
      }
    }
  }
  align_sums
}

fn get_match_probs(
  align_sums: &AlignSums,
  seq_len_pair: &(usize, usize),
  align_scores: &AlignScores,
) -> ProbMat {
  let mut match_probs = vec![vec![0.; seq_len_pair.1]; seq_len_pair.0];
  let mut global_sum = align_sums.forward_sums_match[seq_len_pair.0 - 2][seq_len_pair.1 - 2];
  logsumexp(
    &mut global_sum,
    align_sums.forward_sums_insert[seq_len_pair.0 - 2][seq_len_pair.1 - 2],
  );
  logsumexp(
    &mut global_sum,
    align_sums.forward_sums_del[seq_len_pair.0 - 2][seq_len_pair.1 - 2],
  );
  for (i, x) in match_probs.iter_mut().enumerate() {
    if i == 0 || i == seq_len_pair.0 - 1 {
      continue;
    }
    for (j, x) in x.iter_mut().enumerate() {
      if j == 0 || j == seq_len_pair.1 - 1 {
        continue;
      }
      let mut sum = NEG_INFINITY;
      let forward_sum = align_sums.forward_sums_match[i][j];
      let ends_sum = (i + 1, j + 1) == (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
      let term = if ends_sum {
        0.
      } else {
        align_scores.match2match_score
      } + align_sums.backward_sums_match[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let term = align_scores.match2insert_score + align_sums.backward_sums_insert[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let term = align_scores.match2insert_score + align_sums.backward_sums_del[i + 1][j + 1];
      logsumexp(&mut sum, term);
      let match_prob = expf(forward_sum + sum - global_sum);
      *x = match_prob;
    }
  }
  match_probs
}
