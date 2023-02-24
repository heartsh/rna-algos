use utils::*;

pub struct FoldSums<T: Hash> {
  pub sums_external: SumMat,
  pub sums_rightmost_basepairs_external: SumMat,
  pub sums_rightmost_basepairs_multibranch: SumMat,
  pub sums_close: SparseSumMat<T>,
  pub sums_accessible: SparseSumMat<T>,
  pub sums_multibranch: SumMat,
  pub sums_1ormore_basepairs: SumMat,
}

#[derive(Clone)]
pub struct FoldScores<T: Hash> {
  pub hairpin_scores: SparseScoreMat<T>,
  pub twoloop_scores: ScoreMat4d<T>,
  pub multibranch_close_scores: SparseScoreMat<T>,
  pub accessible_scores: SparseScoreMat<T>,
}

pub type ScoreMat4d<T> = HashMap<PosQuad<T>, Score>;
pub type SparseScoreMat<T> = HashMap<PosPair<T>, Score>;

impl FoldScoreSets {
  pub fn new(init_val: Score) -> FoldScoreSets {
    let init_vals = [init_val; NUM_BASES];
    let mat_2d = [init_vals; NUM_BASES];
    let mat_3d = [mat_2d; NUM_BASES];
    let mat_4d = [mat_3d; NUM_BASES];
    FoldScoreSets {
      // The CONTRAfold model.
      hairpin_scores_len: [init_val; MAX_LOOP_LEN + 1],
      bulge_scores_len: [init_val; MAX_LOOP_LEN],
      interior_scores_len: [init_val; MAX_LOOP_LEN - 1],
      interior_scores_symmetric: [init_val; MAX_INTERIOR_SYMMETRIC],
      interior_scores_asymmetric: [init_val; MAX_INTERIOR_ASYMMETRIC],
      stack_scores: mat_4d,
      terminal_mismatch_scores: mat_4d,
      dangling_scores_left: mat_3d,
      dangling_scores_right: mat_3d,
      helix_close_scores: mat_2d,
      basepair_scores: mat_2d,
      interior_scores_explicit: [[init_val; MAX_INTERIOR_EXPLICIT]; MAX_INTERIOR_EXPLICIT],
      bulge_scores_0x1: [init_val; NUM_BASES],
      interior_scores_1x1: [[init_val; NUM_BASES]; NUM_BASES],
      multibranch_score_base: init_val,
      multibranch_score_basepair: init_val,
      multibranch_score_unpair: init_val,
      external_score_basepair: init_val,
      external_score_unpair: init_val,
      // The cumulative parameters of the CONTRAfold model.
      hairpin_scores_len_cumulative: [init_val; MAX_LOOP_LEN + 1],
      bulge_scores_len_cumulative: [init_val; MAX_LOOP_LEN],
      interior_scores_len_cumulative: [init_val; MAX_LOOP_LEN - 1],
      interior_scores_symmetric_cumulative: [init_val; MAX_INTERIOR_SYMMETRIC],
      interior_scores_asymmetric_cumulative: [init_val; MAX_INTERIOR_ASYMMETRIC],
    }
  }

  pub fn accumulate(&mut self) {
    let mut sum = 0.;
    for i in 0..self.hairpin_scores_len_cumulative.len() {
      sum += self.hairpin_scores_len[i];
      self.hairpin_scores_len_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0..self.bulge_scores_len_cumulative.len() {
      sum += self.bulge_scores_len[i];
      self.bulge_scores_len_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0..self.interior_scores_len_cumulative.len() {
      sum += self.interior_scores_len[i];
      self.interior_scores_len_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0..self.interior_scores_symmetric_cumulative.len() {
      sum += self.interior_scores_symmetric[i];
      self.interior_scores_symmetric_cumulative[i] = sum;
    }
    let mut sum = 0.;
    for i in 0..self.interior_scores_asymmetric_cumulative.len() {
      sum += self.interior_scores_asymmetric[i];
      self.interior_scores_asymmetric_cumulative[i] = sum;
    }
  }

  pub fn transfer(&mut self) {
    for (v, &w) in self
      .hairpin_scores_len
      .iter_mut()
      .zip(HAIRPIN_SCORES_LEN_ATLEAST.iter())
    {
      *v = w;
    }
    for (v, &w) in self
      .bulge_scores_len
      .iter_mut()
      .zip(BULGE_SCORES_LEN_ATLEAST.iter())
    {
      *v = w;
    }
    for (v, &w) in self
      .interior_scores_len
      .iter_mut()
      .zip(INTERIOR_SCORES_LEN_ATLEAST.iter())
    {
      *v = w;
    }
    for (v, &w) in self
      .interior_scores_symmetric
      .iter_mut()
      .zip(INTERIOR_SCORES_SYMMETRIC_ATLEAST.iter())
    {
      *v = w;
    }
    for (v, &w) in self
      .interior_scores_asymmetric
      .iter_mut()
      .zip(INTERIOR_SCORES_ASYMMETRIC_ATLEAST.iter())
    {
      *v = w;
    }
    for (i, x) in STACK_SCORES_CONTRA.iter().enumerate() {
      for (j, x) in x.iter().enumerate() {
        if !has_canonical_basepair(&(i, j)) {
          continue;
        }
        for (k, x) in x.iter().enumerate() {
          for (l, &x) in x.iter().enumerate() {
            if !has_canonical_basepair(&(k, l)) {
              continue;
            }
            self.stack_scores[i][j][k][l] = x;
          }
        }
      }
    }
    for (i, x) in TERMINAL_MISMATCH_SCORES_CONTRA.iter().enumerate() {
      for (j, x) in x.iter().enumerate() {
        if !has_canonical_basepair(&(i, j)) {
          continue;
        }
        for (k, x) in x.iter().enumerate() {
          for (l, &x) in x.iter().enumerate() {
            self.terminal_mismatch_scores[i][j][k][l] = x;
          }
        }
      }
    }
    for (i, x) in DANGLING_SCORES_LEFT.iter().enumerate() {
      for (j, x) in x.iter().enumerate() {
        if !has_canonical_basepair(&(i, j)) {
          continue;
        }
        for (k, &x) in x.iter().enumerate() {
          self.dangling_scores_left[i][j][k] = x;
        }
      }
    }
    for (i, x) in DANGLING_SCORES_RIGHT.iter().enumerate() {
      for (j, x) in x.iter().enumerate() {
        if !has_canonical_basepair(&(i, j)) {
          continue;
        }
        for (k, &x) in x.iter().enumerate() {
          self.dangling_scores_right[i][j][k] = x;
        }
      }
    }
    for (i, x) in HELIX_CLOSE_SCORES.iter().enumerate() {
      for (j, &x) in x.iter().enumerate() {
        if !has_canonical_basepair(&(i, j)) {
          continue;
        }
        self.helix_close_scores[i][j] = x;
      }
    }
    for (i, x) in BASEPAIR_SCORES.iter().enumerate() {
      for (j, &x) in x.iter().enumerate() {
        if !has_canonical_basepair(&(i, j)) {
          continue;
        }
        self.basepair_scores[i][j] = x;
      }
    }
    for (i, x) in INTERIOR_SCORES_EXPLICIT.iter().enumerate() {
      for (j, &x) in x.iter().enumerate() {
        self.interior_scores_explicit[i][j] = x;
      }
    }
    for (x, &y) in self
      .bulge_scores_0x1
      .iter_mut()
      .zip(BULGE_SCORES_0X1.iter())
    {
      *x = y;
    }
    for (i, x) in INTERIOR_SCORES_1X1_CONTRA.iter().enumerate() {
      for (j, &x) in x.iter().enumerate() {
        self.interior_scores_1x1[i][j] = x;
      }
    }
    self.multibranch_score_base = MULTIBRANCH_SCORE_BASE;
    self.multibranch_score_basepair = MULTIBRANCH_SCORE_BASEPAIR;
    self.multibranch_score_unpair = MULTIBRANCH_SCORE_UNPAIR;
    self.external_score_basepair = EXTERNAL_SCORE_BASEPAIR;
    self.external_score_unpair = EXTERNAL_SCORE_UNPAIR;
    self.accumulate();
  }
}

impl<T: Hash> FoldSums<T> {
  pub fn new(seq_len: usize) -> FoldSums<T> {
    let neg_infs = vec![vec![NEG_INFINITY; seq_len]; seq_len];
    FoldSums {
      sums_external: vec![vec![0.; seq_len]; seq_len],
      sums_rightmost_basepairs_external: neg_infs.clone(),
      sums_rightmost_basepairs_multibranch: neg_infs.clone(),
      sums_close: SparseSumMat::<T>::default(),
      sums_accessible: SparseSumMat::<T>::default(),
      sums_multibranch: neg_infs.clone(),
      sums_1ormore_basepairs: neg_infs,
    }
  }
}

impl<T: HashIndex> Default for FoldScores<T> {
  fn default() -> Self {
    Self::new()
  }
}

impl<T: HashIndex> FoldScores<T> {
  pub fn new() -> FoldScores<T> {
    let scores = SparseScoreMat::<T>::default();
    let scores_4d = ScoreMat4d::<T>::default();
    FoldScores {
      hairpin_scores: scores.clone(),
      twoloop_scores: scores_4d,
      multibranch_close_scores: scores.clone(),
      accessible_scores: scores,
    }
  }
}

pub fn mccaskill_algo<T>(
  seq: SeqSlice,
  uses_contra_model: bool,
  allows_short_hairpins: bool,
  fold_score_sets: &FoldScoreSets,
) -> (SparseProbMat<T>, FoldScores<T>)
where
  T: HashIndex,
{
  let seq_len = seq.len();
  let mut fold_scores = FoldScores::<T>::new();
  let fold_sums = if uses_contra_model {
    get_fold_sums_contra::<T>(
      seq,
      &mut fold_scores,
      allows_short_hairpins,
      fold_score_sets,
    )
  } else {
    get_fold_sums::<T>(seq, &mut fold_scores)
  };
  let basepair_probs = if uses_contra_model {
    get_basepair_probs_contra::<T>(
      &fold_sums,
      seq_len,
      &fold_scores,
      allows_short_hairpins,
      fold_score_sets,
    )
  } else {
    get_basepair_probs::<T>(&fold_sums, seq_len, &fold_scores)
  };
  (basepair_probs, fold_scores)
}

pub fn get_fold_sums<T>(seq: SeqSlice, fold_scores: &mut FoldScores<T>) -> FoldSums<T>
where
  T: HashIndex,
{
  let seq_len = seq.len();
  let uses_sentinel_bases = false;
  let mut fold_sums = FoldSums::<T>::new(seq_len);
  let seq_len = T::from_usize(seq.len()).unwrap();
  for subseq_len in range_inclusive(T::from_usize(MIN_SPAN_HAIRPIN_CLOSE).unwrap(), seq_len) {
    for i in range_inclusive(T::zero(), seq_len - subseq_len) {
      let j = i + subseq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let pos_pair_close = (i, j);
      let long_pos_pair_close = (long_i, long_j);
      let basepair_close = (seq[long_i], seq[long_j]);
      let mut sum = NEG_INFINITY;
      if long_pos_pair_close.1 - long_pos_pair_close.0 + 1 >= MIN_SPAN_HAIRPIN_CLOSE
        && has_canonical_basepair(&basepair_close)
      {
        let hairpin_score = get_hairpin_score(seq, &long_pos_pair_close);
        fold_scores
          .hairpin_scores
          .insert(pos_pair_close, hairpin_score);
        logsumexp(&mut sum, hairpin_score);
        for k in range(i + T::one(), j - T::one()) {
          let long_k = k.to_usize().unwrap();
          if long_k - long_i - 1 > MAX_2LOOP_LEN {
            break;
          }
          for l in range(k + T::one(), j).rev() {
            let long_l = l.to_usize().unwrap();
            if long_j - long_l - 1 + long_k - long_i - 1 > MAX_2LOOP_LEN {
              break;
            }
            let pos_pair_accessible = (k, l);
            let long_pos_pair_accessible = (long_k, long_l);
            if let Some(&x) = fold_sums.sums_close.get(&pos_pair_accessible) {
              let y = get_2loop_score(seq, &long_pos_pair_close, &long_pos_pair_accessible);
              fold_scores.twoloop_scores.insert((i, j, k, l), y);
              let y = x + y;
              logsumexp(&mut sum, y);
            }
          }
        }
        let multibranch_close_score = get_multibranch_close_score(seq, &long_pos_pair_close);
        logsumexp(
          &mut sum,
          fold_sums.sums_multibranch[long_i + 1][long_j - 1] + multibranch_close_score,
        );
        let accessible_score = get_accessible_score(seq, &long_pos_pair_close, uses_sentinel_bases);
        if sum > NEG_INFINITY {
          fold_scores
            .multibranch_close_scores
            .insert(pos_pair_close, multibranch_close_score);
          fold_scores
            .accessible_scores
            .insert(pos_pair_close, accessible_score);
          fold_sums.sums_close.insert(pos_pair_close, sum);
          let sum = sum + accessible_score;
          fold_sums.sums_accessible.insert(pos_pair_close, sum);
        }
      }
      sum = NEG_INFINITY;
      for k in range_inclusive(i + T::one(), j) {
        let pos_pair_accessible = (i, k);
        if let Some(&x) = fold_sums.sums_accessible.get(&pos_pair_accessible) {
          logsumexp(&mut sum, x);
        }
      }
      fold_sums.sums_rightmost_basepairs_external[long_i][long_j] = sum;
      sum = 0.;
      for k in long_i..long_j {
        let x = fold_sums.sums_rightmost_basepairs_external[k][long_j];
        let y = if long_i == 0 && k == 0 {
          0.
        } else {
          fold_sums.sums_external[long_i][k - 1]
        };
        let y = x + y;
        logsumexp(&mut sum, y);
      }
      fold_sums.sums_external[long_i][long_j] = sum;
      sum = fold_sums.sums_rightmost_basepairs_external[long_i][long_j] + COEFF_NUM_BRANCHES;
      let mut sum2 = NEG_INFINITY;
      for k in long_i + 1..long_j {
        let x = fold_sums.sums_rightmost_basepairs_external[k][long_j] + COEFF_NUM_BRANCHES;
        logsumexp(&mut sum, x);
        let y = fold_sums.sums_1ormore_basepairs[long_i][k - 1] + x;
        logsumexp(&mut sum2, y);
      }
      fold_sums.sums_multibranch[long_i][long_j] = sum2;
      logsumexp(&mut sum, sum2);
      fold_sums.sums_1ormore_basepairs[long_i][long_j] = sum;
    }
  }
  fold_sums
}

pub fn get_fold_sums_contra<T>(
  seq: SeqSlice,
  fold_scores: &mut FoldScores<T>,
  allows_short_hairpins: bool,
  fold_score_sets: &FoldScoreSets,
) -> FoldSums<T>
where
  T: HashIndex,
{
  let seq_len = seq.len();
  let uses_sentinel_bases = false;
  let mut fold_sums = FoldSums::<T>::new(seq_len);
  let short_seq_len = T::from_usize(seq.len()).unwrap();
  for subseq_len in range_inclusive(T::one(), short_seq_len) {
    for i in range_inclusive(T::zero(), short_seq_len - subseq_len) {
      let j = i + subseq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let pos_pair_close = (i, j);
      let long_pos_pair_close = (long_i, long_j);
      let basepair_close = (seq[long_i], seq[long_j]);
      let mut sum = NEG_INFINITY;
      if has_canonical_basepair(&basepair_close)
        && (allows_short_hairpins
          || long_pos_pair_close.1 - long_pos_pair_close.0 + 1 >= MIN_SPAN_HAIRPIN_CLOSE)
      {
        if long_j - long_i - 1 <= MAX_LOOP_LEN {
          let hairpin_score = get_hairpin_score_contra(seq, &long_pos_pair_close, fold_score_sets);
          fold_scores
            .hairpin_scores
            .insert(pos_pair_close, hairpin_score);
          logsumexp(&mut sum, hairpin_score);
        }
        for k in range(i + T::one(), j - T::one()) {
          let long_k = k.to_usize().unwrap();
          if long_k - long_i - 1 > MAX_LOOP_LEN {
            break;
          }
          for l in range(k + T::one(), j).rev() {
            let long_l = l.to_usize().unwrap();
            if long_j - long_l - 1 + long_k - long_i - 1 > MAX_LOOP_LEN {
              break;
            }
            let pos_pair_accessible = (k, l);
            let long_pos_pair_accessible = (long_k, long_l);
            if let Some(&x) = fold_sums.sums_close.get(&pos_pair_accessible) {
              let y = get_2loop_score_contra(
                seq,
                &long_pos_pair_close,
                &long_pos_pair_accessible,
                fold_score_sets,
              );
              fold_scores.twoloop_scores.insert((i, j, k, l), y);
              let y = x + y;
              logsumexp(&mut sum, y);
            }
          }
        }
        let multibranch_close_score = fold_score_sets.multibranch_score_base
          + fold_score_sets.multibranch_score_basepair
          + get_junction_score(
            seq,
            &long_pos_pair_close,
            uses_sentinel_bases,
            fold_score_sets,
          );
        logsumexp(
          &mut sum,
          fold_sums.sums_multibranch[long_i + 1][long_j - 1] + multibranch_close_score,
        );
        let accessible_score = get_junction_score(
          seq,
          &(long_pos_pair_close.1, long_pos_pair_close.0),
          uses_sentinel_bases,
          fold_score_sets,
        ) + fold_score_sets.basepair_scores[basepair_close.0]
          [basepair_close.1];
        if sum > NEG_INFINITY {
          fold_scores
            .multibranch_close_scores
            .insert(pos_pair_close, multibranch_close_score);
          fold_scores
            .accessible_scores
            .insert(pos_pair_close, accessible_score);
          fold_sums.sums_close.insert(pos_pair_close, sum);
          let sum = sum + accessible_score;
          fold_sums.sums_accessible.insert(pos_pair_close, sum);
        }
      }
      sum = NEG_INFINITY;
      let mut sum2 = sum;
      for k in range_inclusive(i + T::one(), j) {
        let pos_pair_accessible = (i, k);
        if let Some(&x) = fold_sums.sums_accessible.get(&pos_pair_accessible) {
          logsumexp(
            &mut sum,
            x + fold_score_sets.external_score_basepair
              + fold_score_sets.external_score_unpair * (j - k).to_f32().unwrap(),
          );
          logsumexp(
            &mut sum2,
            x + fold_score_sets.multibranch_score_basepair
              + fold_score_sets.multibranch_score_unpair * (j - k).to_f32().unwrap(),
          );
        }
      }
      fold_sums.sums_rightmost_basepairs_external[long_i][long_j] = sum;
      fold_sums.sums_rightmost_basepairs_multibranch[long_i][long_j] = sum2;
      sum = fold_score_sets.external_score_unpair * subseq_len.to_f32().unwrap();
      for k in long_i..long_j {
        let x = fold_sums.sums_rightmost_basepairs_external[k][long_j];
        let y = if long_i == 0 && k == 0 {
          0.
        } else {
          fold_sums.sums_external[long_i][k - 1]
        };
        let y = x + y;
        logsumexp(&mut sum, y);
      }
      fold_sums.sums_external[long_i][long_j] = sum;
      sum = fold_sums.sums_rightmost_basepairs_multibranch[long_i][long_j];
      sum2 = NEG_INFINITY;
      for k in long_i + 1..long_j {
        let x = fold_sums.sums_rightmost_basepairs_multibranch[k][long_j];
        logsumexp(
          &mut sum,
          x + fold_score_sets.multibranch_score_unpair * (k - long_i) as Score,
        );
        let x = fold_sums.sums_1ormore_basepairs[long_i][k - 1] + x;
        logsumexp(&mut sum2, x);
      }
      fold_sums.sums_multibranch[long_i][long_j] = sum2;
      logsumexp(&mut sum, sum2);
      fold_sums.sums_1ormore_basepairs[long_i][long_j] = sum;
    }
  }
  fold_sums
}

pub fn get_basepair_probs<T>(
  fold_sums: &FoldSums<T>,
  seq_len: usize,
  fold_scores: &FoldScores<T>,
) -> SparseProbMat<T>
where
  T: HashIndex,
{
  let global_sum = fold_sums.sums_external[0][seq_len - 1];
  let mut basepair_probs = SparseProbMat::<T>::default();
  let mut probs_multibranch = vec![vec![NEG_INFINITY; seq_len]; seq_len];
  let mut probs_multibranch2 = probs_multibranch.clone();
  let short_seq_len = T::from_usize(seq_len).unwrap();
  for subseq_len in range_inclusive(
    T::from_usize(MIN_SPAN_HAIRPIN_CLOSE).unwrap(),
    short_seq_len,
  )
  .rev()
  {
    for i in range_inclusive(T::zero(), short_seq_len - subseq_len) {
      let j = i + subseq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let mut sum = NEG_INFINITY;
      let mut sum2 = sum;
      for k in range(j + T::one(), short_seq_len) {
        let long_k = k.to_usize().unwrap();
        let pos_pair_close = (i, k);
        if let Some(&x) = fold_sums.sums_close.get(&pos_pair_close) {
          let basepair_prob = basepair_probs[&pos_pair_close];
          let multibranch_close_score = fold_scores.multibranch_close_scores[&pos_pair_close];
          let x = basepair_prob + multibranch_close_score - x;
          logsumexp(
            &mut sum,
            x + fold_sums.sums_1ormore_basepairs[long_j + 1][long_k - 1],
          );
          logsumexp(&mut sum2, x);
        }
      }
      probs_multibranch[long_i][long_j] = sum;
      probs_multibranch2[long_i][long_j] = sum2;
      let pos_pair_accessible = (i, j);
      if let Some(&sum_close) = fold_sums.sums_close.get(&pos_pair_accessible) {
        let sum_accessible = fold_sums.sums_accessible[&pos_pair_accessible];
        let sum_pair = (
          if pos_pair_accessible.0 < T::one() {
            0.
          } else {
            fold_sums.sums_external[0][long_i - 1]
          },
          if pos_pair_accessible.1 > short_seq_len - T::from_usize(2).unwrap() {
            0.
          } else {
            fold_sums.sums_external[long_j + 1][seq_len - 1]
          },
        );
        let mut sum = sum_pair.0 + sum_accessible + sum_pair.1 - global_sum;
        for k in range(T::zero(), i).rev() {
          let long_k = k.to_usize().unwrap();
          if long_i - long_k - 1 > MAX_2LOOP_LEN {
            break;
          }
          for l in range(j + T::one(), short_seq_len) {
            let long_l = l.to_usize().unwrap();
            if long_l - long_j - 1 + long_i - long_k - 1 > MAX_2LOOP_LEN {
              break;
            }
            let pos_pair_close = (k, l);
            if let Some(&x) = fold_sums.sums_close.get(&pos_pair_close) {
              logsumexp(
                &mut sum,
                basepair_probs[&pos_pair_close] + sum_close - x
                  + fold_scores.twoloop_scores[&(k, l, i, j)],
              );
            }
          }
        }
        let sum_accessible = sum_accessible + COEFF_NUM_BRANCHES;
        for k in 0..long_i {
          let x = fold_sums.sums_1ormore_basepairs[k + 1][long_i - 1];
          logsumexp(&mut sum, sum_accessible + probs_multibranch2[k][long_j] + x);
          let y = probs_multibranch[k][long_j];
          logsumexp(&mut sum, sum_accessible + y);
          logsumexp(&mut sum, sum_accessible + x + y);
        }
        if sum > NEG_INFINITY {
          basepair_probs.insert(pos_pair_accessible, sum);
        }
      }
    }
  }
  basepair_probs = basepair_probs.iter().map(|(x, &y)| (*x, expf(y))).collect();
  basepair_probs
}

pub fn get_basepair_probs_contra<T>(
  fold_sums: &FoldSums<T>,
  seq_len: usize,
  fold_scores: &FoldScores<T>,
  allows_short_hairpins: bool,
  fold_score_sets: &FoldScoreSets,
) -> SparseProbMat<T>
where
  T: HashIndex,
{
  let global_sum = fold_sums.sums_external[0][seq_len - 1];
  let mut basepair_probs = SparseProbMat::<T>::default();
  let mut probs_multibranch = vec![vec![NEG_INFINITY; seq_len]; seq_len];
  let mut probs_multibranch2 = probs_multibranch.clone();
  let short_seq_len = T::from_usize(seq_len).unwrap();
  for subseq_len in range_inclusive(
    T::from_usize(if allows_short_hairpins {
      2
    } else {
      MIN_SPAN_HAIRPIN_CLOSE
    })
    .unwrap(),
    short_seq_len,
  )
  .rev()
  {
    for i in range_inclusive(T::zero(), short_seq_len - subseq_len) {
      let j = i + subseq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let mut sum = NEG_INFINITY;
      let mut sum2 = sum;
      for k in range(j + T::one(), short_seq_len) {
        let long_k = k.to_usize().unwrap();
        let pos_pair_close = (i, k);
        if let Some(&x) = fold_sums.sums_close.get(&pos_pair_close) {
          let basepair_prob = basepair_probs[&pos_pair_close];
          let multibranch_close_score = fold_scores.multibranch_close_scores[&pos_pair_close];
          let x = basepair_prob + multibranch_close_score - x;
          logsumexp(
            &mut sum,
            x + fold_sums.sums_1ormore_basepairs[long_j + 1][long_k - 1],
          );
          logsumexp(
            &mut sum2,
            x + fold_score_sets.multibranch_score_unpair * (k - j - T::one()).to_f32().unwrap(),
          );
        }
      }
      probs_multibranch[long_i][long_j] = sum;
      probs_multibranch2[long_i][long_j] = sum2;
      let pos_pair_accessible = (i, j);
      if let Some(&sum_close) = fold_sums.sums_close.get(&pos_pair_accessible) {
        let sum_pair = (
          if pos_pair_accessible.0 < T::one() {
            0.
          } else {
            fold_sums.sums_external[0][long_i - 1]
          },
          if pos_pair_accessible.1 > short_seq_len - T::from_usize(2).unwrap() {
            0.
          } else {
            fold_sums.sums_external[long_j + 1][seq_len - 1]
          },
        );
        let mut sum = sum_pair.0
          + sum_pair.1
          + fold_sums.sums_accessible[&pos_pair_accessible]
          + fold_score_sets.external_score_basepair
          - global_sum;
        for k in range(T::zero(), i).rev() {
          let long_k = k.to_usize().unwrap();
          if long_i - long_k - 1 > MAX_LOOP_LEN {
            break;
          }
          for l in range(j + T::one(), short_seq_len) {
            let long_l = l.to_usize().unwrap();
            if long_l - long_j - 1 + long_i - long_k - 1 > MAX_LOOP_LEN {
              break;
            }
            let pos_pair_close = (k, l);
            if let Some(&x) = fold_sums.sums_close.get(&pos_pair_close) {
              logsumexp(
                &mut sum,
                basepair_probs[&pos_pair_close] + sum_close - x
                  + fold_scores.twoloop_scores[&(k, l, i, j)],
              );
            }
          }
        }
        let sum_accessible = fold_sums.sums_accessible[&pos_pair_accessible]
          + fold_score_sets.multibranch_score_basepair;
        for k in 0..long_i {
          let x = fold_sums.sums_1ormore_basepairs[k + 1][long_i - 1];
          logsumexp(&mut sum, sum_accessible + probs_multibranch2[k][long_j] + x);
          let y = probs_multibranch[k][long_j];
          logsumexp(
            &mut sum,
            sum_accessible
              + y
              + fold_score_sets.multibranch_score_unpair * (long_i - k - 1) as Score,
          );
          logsumexp(&mut sum, sum_accessible + x + y);
        }
        if sum > NEG_INFINITY {
          basepair_probs.insert(pos_pair_accessible, sum);
        }
      }
    }
  }
  basepair_probs = basepair_probs.iter().map(|(x, &y)| (*x, expf(y))).collect();
  basepair_probs
}
