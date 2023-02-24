use utils::*;

#[derive(Clone)]
pub struct CentroidFold<T> {
  pub basepair_pos_pairs: PosPairs<T>,
  pub expect_accuracy: Score,
}
pub type Poss<T> = Vec<T>;

impl<T> Default for CentroidFold<T> {
  fn default() -> Self {
    Self::new()
  }
}

impl<T> CentroidFold<T> {
  pub fn new() -> CentroidFold<T> {
    CentroidFold {
      basepair_pos_pairs: PosPairs::<T>::new(),
      expect_accuracy: 0.,
    }
  }
}

pub fn centroid_fold<T: HashIndex>(
  basepair_probs: &SparseProbMat<T>,
  seq_len: usize,
  centroid_threshold: Prob,
) -> CentroidFold<T>
where
  T: HashIndex,
{
  let mut max_expect_accuracies = vec![vec![0.; seq_len]; seq_len];
  let seq_len = T::from_usize(seq_len).unwrap();
  for subseq_len in range_inclusive(T::one(), seq_len) {
    for i in range_inclusive(T::zero(), seq_len - subseq_len) {
      let j = i + subseq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      if i == j {
        continue;
      }
      let mut max_expect_accuracy = max_expect_accuracies[long_i + 1][long_j];
      let expect_accuracy = max_expect_accuracies[long_i][long_j - 1];
      if expect_accuracy > max_expect_accuracy {
        max_expect_accuracy = expect_accuracy;
      }
      let pos_pair = (i, j);
      if let Some(&x) = basepair_probs.get(&pos_pair) {
        let expect_accuracy = max_expect_accuracies[long_i + 1][long_j - 1] + centroid_threshold * x - 1.;
        if expect_accuracy > max_expect_accuracy {
          max_expect_accuracy = expect_accuracy;
        }
      }
      for k in long_i + 1..long_j {
        let expect_accuracy = max_expect_accuracies[long_i][k] + max_expect_accuracies[k + 1][long_j];
        if expect_accuracy > max_expect_accuracy {
          max_expect_accuracy = expect_accuracy;
        }
      }
      max_expect_accuracies[long_i][long_j] = max_expect_accuracy;
    }
  }
  let mut centroid_fold = CentroidFold::<T>::new();
  let mut pos_pair_stack = vec![(T::zero(), seq_len - T::one())];
  while !pos_pair_stack.is_empty() {
    let pos_pair = pos_pair_stack.pop().unwrap();
    let (i, j) = pos_pair;
    if j <= i {
      continue;
    }
    let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
    let max_expect_accuracy = max_expect_accuracies[long_i][long_j];
    if max_expect_accuracy == 0. {
      continue;
    }
    if max_expect_accuracy == max_expect_accuracies[long_i + 1][long_j] {
      pos_pair_stack.push((i + T::one(), j));
    } else if max_expect_accuracy == max_expect_accuracies[long_i][long_j - 1] {
      pos_pair_stack.push((i, j - T::one()));
    } else if basepair_probs.contains_key(&pos_pair)
      && max_expect_accuracy == max_expect_accuracies[long_i + 1][long_j - 1] + centroid_threshold * basepair_probs[&pos_pair] - 1.
    {
      pos_pair_stack.push((i + T::one(), j - T::one()));
      centroid_fold.basepair_pos_pairs.push(pos_pair);
    } else {
      for k in range(i + T::one(), j) {
        let long_k = k.to_usize().unwrap();
        if max_expect_accuracy == max_expect_accuracies[long_i][long_k] + max_expect_accuracies[long_k + 1][long_j] {
          pos_pair_stack.push((i, k));
          pos_pair_stack.push((k + T::one(), j));
          break;
        }
      }
    }
  }
  centroid_fold.expect_accuracy = max_expect_accuracies[0][seq_len.to_usize().unwrap() - 1];
  centroid_fold
}
