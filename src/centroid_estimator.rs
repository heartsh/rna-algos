use utils::*;

#[derive(Clone)]
pub struct MeaSs<T> {
  pub bp_pos_pairs: PosPairs<T>,
  pub ea: Mea,
}
pub type Meas = Vec<Mea>;
pub type MeaMat<T> = HashMap<PosPair<T>, Mea>;
pub type Poss<T> = Vec<T>;
pub type PosSeqsWithPoss<T> = HashMap<T, Poss<T>>;
pub type PosPairSeqsWithPosPairs<T> = HashMap<PosPair<T>, PosPairs<T>>;

impl<T> MeaSs<T> {
  pub fn new() -> MeaSs<T> {
    MeaSs {
      bp_pos_pairs: PosPairs::<T>::new(),
      ea: 0.,
    }
  }
}

pub fn centroid_estimator<T: HashIndex>(
  bpp_mat: &SparseProbMat<T>,
  seq_len: usize,
  gamma: Prob,
) -> MeaSs<T>
where
  T: HashIndex,
{
  let gamma_plus_1 = gamma + 1.;
  let mut mea_mat = vec![vec![0.; seq_len]; seq_len];
  let seq_len = T::from_usize(seq_len).unwrap();
  for sub_seq_len in range_inclusive(T::one(), seq_len) {
    for i in range_inclusive(T::zero(), seq_len - sub_seq_len) {
      let j = i + sub_seq_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      if i == j {
        continue;
      }
      let mut mea = mea_mat[long_i + 1][long_j];
      let ea = mea_mat[long_i][long_j - 1];
      if ea > mea {
        mea = ea;
      }
      let pos_pair = (i, j);
      match bpp_mat.get(&pos_pair) {
        Some(&bpp) => {
          let ea = mea_mat[long_i + 1][long_j - 1] + gamma_plus_1 * bpp - 1.;
          if ea > mea {
            mea = ea;
          }
        }
        None => {}
      }
      for k in long_i + 1..long_j {
        let ea = mea_mat[long_i][k] + mea_mat[k + 1][long_j];
        if ea > mea {
          mea = ea;
        }
      }
      mea_mat[long_i][long_j] = mea;
    }
  }
  let mut mea_ss = MeaSs::<T>::new();
  let mut pos_pair_stack = vec![(T::zero(), seq_len - T::one())];
  while pos_pair_stack.len() > 0 {
    let pos_pair = pos_pair_stack.pop().unwrap();
    let (i, j) = pos_pair;
    if j <= i {
      continue;
    }
    let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
    let mea = mea_mat[long_i][long_j];
    if mea == 0. {
      continue;
    }
    if mea == mea_mat[long_i + 1][long_j] {
      pos_pair_stack.push((i + T::one(), j));
    } else if mea == mea_mat[long_i][long_j - 1] {
      pos_pair_stack.push((i, j - T::one()));
    } else if bpp_mat.contains_key(&pos_pair)
      && mea == mea_mat[long_i + 1][long_j - 1] + gamma_plus_1 * bpp_mat[&pos_pair] - 1.
    {
      pos_pair_stack.push((i + T::one(), j - T::one()));
      mea_ss.bp_pos_pairs.push(pos_pair);
    } else {
      for k in range(i + T::one(), j) {
        let long_k = k.to_usize().unwrap();
        if mea == mea_mat[long_i][long_k] + mea_mat[long_k + 1][long_j] {
          pos_pair_stack.push((i, k));
          pos_pair_stack.push((k + T::one(), j));
          break;
        }
      }
    }
  }
  mea_ss.ea = mea_mat[0][seq_len.to_usize().unwrap() - 1];
  mea_ss
}
