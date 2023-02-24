pub use bio::io::fasta::Reader;
pub use hashbrown::HashMap;
pub use itertools::multizip;
pub use num::{
  range, range_inclusive, Bounded, FromPrimitive, Integer, One, PrimInt, ToPrimitive, Unsigned,
  Zero,
};
pub use rna_ss_params::compiled_scores_contra::*;
pub use rna_ss_params::compiled_scores_turner::*;
pub use rna_ss_params::utils::*;
pub use scoped_threadpool::Pool;
pub use std::cmp::{max, min};
pub use std::env;
pub use std::f32::NEG_INFINITY;
pub use std::fmt::Display;
pub use std::fs::create_dir;
pub use std::fs::File;
pub use std::hash::Hash;
pub use std::io::prelude::*;
pub use std::io::{BufReader, BufWriter};
pub use std::path::Path;
pub use std::str::from_utf8_unchecked;

pub trait HashIndex:
  Unsigned
  + PrimInt
  + Hash
  + FromPrimitive
  + ToPrimitive
  + Clone
  + Integer
  + Eq
  + One
  + Ord
  + Display
  + Sync
  + Send
{
}
impl<T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + One + Ord + Display + Sync + Send>
  HashIndex for T
{
}

pub type PosPair<T> = (T, T);
pub type PosQuad<T> = (T, T, T, T);
type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
#[derive(Clone)]
pub struct FastaRecord {
  pub fasta_id: FastaId,
  pub seq: Seq,
}
#[derive(Debug)]
pub struct Align<T> {
  pub cols: Cols,
  pub pos_map_sets: PosMapSets<T>,
}
pub type PosMaps<T> = Vec<T>;
pub type PosMapSets<T> = Vec<PosMaps<T>>;
pub type FastaRecords = Vec<FastaRecord>;
pub type SeqSlice<'a> = &'a [Base];
pub type NumThreads = u32;
pub type SparseProbMat<T> = HashMap<PosPair<T>, Prob>;
pub type SeqId = String;
pub type SeqIds = Vec<SeqId>;
pub type Col = Vec<Base>;
pub type Cols = Vec<Col>;
pub type Sum4dMat<T> = HashMap<PosQuad<T>, Sum>;
pub type SparseSumMat<T> = HashMap<PosPair<T>, Sum>;
pub type PosPairs<T> = Vec<PosPair<T>>;
pub type FoldChar = u8;
pub type FoldStr = Vec<FoldChar>;
pub type Seq = Vec<Base>;
pub type SeqPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type Prob = f32;
pub type Sum = Prob;
pub type Probs = Vec<Prob>;
pub type ProbMat = Vec<Probs>;
pub type Sums = Vec<Sum>;
pub type SumMat = Vec<Sums>;
pub type Pos = usize;
pub type Base = usize;
pub type MatchScores = [[Score; NUM_BASES]; NUM_BASES];
pub type InsertScores = [Score; NUM_BASES];
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type ProbMatsHashedIds = HashMap<RnaIdPair, ProbMat>;

#[derive(Clone, Debug)]
pub struct FoldScoreSets {
  // The CONTRAfold model.
  pub hairpin_scores_len: HairpinScoresLen,
  pub bulge_scores_len: BulgeScoresLen,
  pub interior_scores_len: InteriorScoresLen,
  pub interior_scores_symmetric: InteriorScoresSymmetric,
  pub interior_scores_asymmetric: InteriorScoresAsymmetric,
  pub stack_scores: StackScores,
  pub terminal_mismatch_scores: TerminalMismatchScores,
  pub dangling_scores_left: DanglingScores,
  pub dangling_scores_right: DanglingScores,
  pub helix_close_scores: HelixCloseScores,
  pub basepair_scores: BasepairScores,
  pub interior_scores_explicit: InteriorScoresExplicit,
  pub bulge_scores_0x1: BulgeScores0x1,
  pub interior_scores_1x1: InteriorScores1x1Contra,
  pub multibranch_score_base: Score,
  pub multibranch_score_basepair: Score,
  pub multibranch_score_unpair: Score,
  pub external_score_basepair: Score,
  pub external_score_unpair: Score,
  // The cumulative parameters of the CONTRAfold model.
  pub hairpin_scores_len_cumulative: HairpinScoresLen,
  pub bulge_scores_len_cumulative: BulgeScoresLen,
  pub interior_scores_len_cumulative: InteriorScoresLen,
  pub interior_scores_symmetric_cumulative: InteriorScoresSymmetric,
  pub interior_scores_asymmetric_cumulative: InteriorScoresAsymmetric,
}

pub const LOGSUMEXP_THRESHOLD_UPPER: Score = 11.862_479;
pub const PSEUDO_BASE: Base = U + 1 as Base;
pub const UNPAIR: FoldChar = b'.';
pub const BASEPAIR_LEFT: FoldChar = b'(';
pub const BASEPAIR_RIGHT: FoldChar = b')';
pub const EXAMPLE_FASTA_FILE_PATH: &str = "assets/sampled_trnas.fa";
pub const EPSILON: Prob = 0.00_1;
pub const PROB_BOUND_LOWER: Prob = -EPSILON;
pub const PROB_BOUND_UPPER: Prob = 1. + EPSILON;

impl FastaRecord {
  pub fn origin() -> FastaRecord {
    FastaRecord {
      fasta_id: FastaId::new(),
      seq: Seq::new(),
    }
  }

  pub fn new(fasta_id: FastaId, seq: Seq) -> FastaRecord {
    FastaRecord {
      fasta_id: fasta_id,
      seq: seq,
    }
  }
}

impl<T> Default for Align<T> {
  fn default() -> Self {
    Self::new()
  }
}

impl<T> Align<T> {
  pub fn new() -> Align<T> {
    Align {
      cols: Cols::new(),
      pos_map_sets: PosMapSets::<T>::new(),
    }
  }
}

pub fn has_canonical_basepair(x: &Basepair) -> bool {
  matches!(*x, AU | CG | GC | GU | UA | UG)
}

pub fn get_hairpin_score(seq: SeqSlice, pos_pair_close: &(usize, usize)) -> Score {
  let hairpin = &seq[pos_pair_close.0..pos_pair_close.1 + 1];
  let special_hairpin_score = get_special_hairpin_score(hairpin);
  if special_hairpin_score > NEG_INFINITY {
    special_hairpin_score
  } else {
    let hairpin_len = pos_pair_close.1 - pos_pair_close.0 - 1;
    let basepair_close = (seq[pos_pair_close.0], seq[pos_pair_close.1]);
    let hairpin_score = if hairpin_len == MIN_HAIRPIN_LEN {
      HAIRPIN_SCORES_INIT[hairpin_len]
    } else {
      let terminal_mismatch = (seq[pos_pair_close.0 + 1], seq[pos_pair_close.1 - 1]);
      let hairpin_score_init = if hairpin_len <= MAX_HAIRPIN_LEN_EXTRAPOLATION {
        HAIRPIN_SCORES_INIT[hairpin_len]
      } else {
        HAIRPIN_SCORES_INIT[MIN_HAIRPIN_LEN_EXTRAPOLATION - 1]
          + COEFF_HAIRPIN_LEN_EXTRAPOLATION
            * (hairpin_len as Score / (MIN_HAIRPIN_LEN_EXTRAPOLATION - 1) as Score).ln()
      };
      hairpin_score_init
        + TERMINAL_MISMATCH_SCORES_HAIRPIN[basepair_close.0][basepair_close.1][terminal_mismatch.0]
          [terminal_mismatch.1]
    };
    hairpin_score
      + if matches_augu(&basepair_close) {
        HELIX_AUGU_END_PENALTY
      } else {
        0.
      }
  }
}

pub fn get_special_hairpin_score(seq: SeqSlice) -> Score {
  for x in HAIRPIN_SCORES_SPECIAL.iter() {
    if x.0 == seq {
      return x.1;
    }
  }
  NEG_INFINITY
}

pub fn get_2loop_score(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
) -> Score {
  if pos_pair_close.0 + 1 == pos_pair_accessible.0 && pos_pair_close.1 - 1 == pos_pair_accessible.1
  {
    get_stack_score(seq, pos_pair_close, pos_pair_accessible)
  } else if pos_pair_close.0 + 1 == pos_pair_accessible.0
    || pos_pair_close.1 - 1 == pos_pair_accessible.1
  {
    get_bulge_score(seq, pos_pair_close, pos_pair_accessible)
  } else {
    get_interior_score(seq, pos_pair_close, pos_pair_accessible)
  }
}

fn get_stack_score(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
) -> Score {
  let basepair_close = (seq[pos_pair_close.0], seq[pos_pair_close.1]);
  let basepair_accessible = (seq[pos_pair_accessible.0], seq[pos_pair_accessible.1]);
  STACK_SCORES[basepair_close.0][basepair_close.1][basepair_accessible.0][basepair_accessible.1]
}

fn get_bulge_score(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
) -> Score {
  let bulge_len =
    pos_pair_accessible.0 - pos_pair_close.0 + pos_pair_close.1 - pos_pair_accessible.1 - 2;
  if bulge_len == 1 {
    BULGE_SCORES_INIT[bulge_len] + get_stack_score(seq, pos_pair_close, pos_pair_accessible)
  } else {
    let basepair_close = (seq[pos_pair_close.0], seq[pos_pair_close.1]);
    let basepair_accessible = (seq[pos_pair_accessible.0], seq[pos_pair_accessible.1]);
    BULGE_SCORES_INIT[bulge_len]
      + if matches_augu(&basepair_close) {
        HELIX_AUGU_END_PENALTY
      } else {
        0.
      }
      + if matches_augu(&basepair_accessible) {
        HELIX_AUGU_END_PENALTY
      } else {
        0.
      }
  }
}

fn get_interior_score(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
) -> Score {
  let basepair_close = (seq[pos_pair_close.0], seq[pos_pair_close.1]);
  let basepair_accessible = (seq[pos_pair_accessible.0], seq[pos_pair_accessible.1]);
  let loop_len_pair = (
    pos_pair_accessible.0 - pos_pair_close.0 - 1,
    pos_pair_close.1 - pos_pair_accessible.1 - 1,
  );
  let interior_len = loop_len_pair.0 + loop_len_pair.1;
  match loop_len_pair {
    (1, 1) => {
      let interior = (seq[pos_pair_close.0 + 1], seq[pos_pair_close.1 - 1]);
      INTERIOR_SCORES_1X1[basepair_close.0][basepair_close.1][interior.0][interior.1]
        [basepair_accessible.0][basepair_accessible.1]
    }
    (1, 2) => {
      let interior = (
        (seq[pos_pair_close.0 + 1], seq[pos_pair_close.1 - 1]),
        seq[pos_pair_close.1 - 2],
      );
      INTERIOR_SCORES_1X2[basepair_close.0][basepair_close.1][(interior.0).0][(interior.0).1]
        [interior.1][basepair_accessible.0][basepair_accessible.1]
    }
    (2, 1) => {
      let interior = (
        (seq[pos_pair_close.1 - 1], seq[pos_pair_close.0 + 2]),
        seq[pos_pair_close.0 + 1],
      );
      let basepair_accessible_inverse = invert_basepair(&basepair_accessible);
      let basepair_close_inverse = invert_basepair(&basepair_close);
      INTERIOR_SCORES_1X2[basepair_accessible_inverse.0][basepair_accessible_inverse.1]
        [(interior.0).0][(interior.0).1][interior.1][basepair_close_inverse.0]
        [basepair_close_inverse.1]
    }
    (2, 2) => {
      let interior = (
        (seq[pos_pair_close.0 + 1], seq[pos_pair_close.1 - 1]),
        (seq[pos_pair_close.0 + 2], seq[pos_pair_close.1 - 2]),
      );
      INTERIOR_SCORES_2X2[basepair_close.0][basepair_close.1][(interior.0).0][(interior.0).1]
        [(interior.1).0][(interior.1).1][basepair_accessible.0][basepair_accessible.1]
    }
    _ => {
      INTERIOR_SCORES_INIT[interior_len]
        + (NINIO_COEFF * get_abs_diff(loop_len_pair.0, loop_len_pair.1) as Score).max(NINIO_MAX)
        + get_interior_mismatch_score(seq, pos_pair_close, pos_pair_accessible, &loop_len_pair)
        + if matches_augu(&basepair_close) {
          HELIX_AUGU_END_PENALTY
        } else {
          0.
        }
        + if matches_augu(&basepair_accessible) {
          HELIX_AUGU_END_PENALTY
        } else {
          0.
        }
    }
  }
}

pub fn invert_basepair(x: &Basepair) -> Basepair {
  (x.1, x.0)
}

pub fn get_abs_diff(x: usize, y: usize) -> usize {
  max(x, y) - min(x, y)
}

fn get_interior_mismatch_score(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
  loop_len_pair: &(usize, usize),
) -> Score {
  let basepair_close = (seq[pos_pair_close.0], seq[pos_pair_close.1]);
  let basepair_accessible = (seq[pos_pair_accessible.1], seq[pos_pair_accessible.0]);
  let terminal_mismatch = (
    (seq[pos_pair_close.0 + 1], seq[pos_pair_close.1 - 1]),
    (
      seq[pos_pair_accessible.1 + 1],
      seq[pos_pair_accessible.0 - 1],
    ),
  );
  match *loop_len_pair {
    (1, _) | (_, 1) => {
      TERMINAL_MISMATCH_SCORES_1XMANY[basepair_close.0][basepair_close.1][(terminal_mismatch.0).0]
        [(terminal_mismatch.0).1]
        + TERMINAL_MISMATCH_SCORES_1XMANY[basepair_accessible.0][basepair_accessible.1]
          [(terminal_mismatch.1).0][(terminal_mismatch.1).1]
    }
    (2, 3) | (3, 2) => {
      TERMINAL_MISMATCH_SCORES_2X3[basepair_close.0][basepair_close.1][(terminal_mismatch.0).0]
        [(terminal_mismatch.0).1]
        + TERMINAL_MISMATCH_SCORES_2X3[basepair_accessible.0][basepair_accessible.1]
          [(terminal_mismatch.1).0][(terminal_mismatch.1).1]
    }
    _ => {
      TERMINAL_MISMATCH_SCORES_INTERIOR[basepair_close.0][basepair_close.1][(terminal_mismatch.0).0]
        [(terminal_mismatch.0).1]
        + TERMINAL_MISMATCH_SCORES_INTERIOR[basepair_accessible.0][basepair_accessible.1]
          [(terminal_mismatch.1).0][(terminal_mismatch.1).1]
    }
  }
}

pub fn get_multibranch_close_score(seq: SeqSlice, pos_pair_close: &(usize, usize)) -> Score {
  let basepair_close = (seq[pos_pair_close.0], seq[pos_pair_close.1]);
  let basepair_close_inverse = invert_basepair(&basepair_close);
  let basepair_stack_inverse =
    invert_basepair(&(seq[pos_pair_close.0 + 1], seq[pos_pair_close.1 - 1]));
  let terminal_mismatch_score = TERMINAL_MISMATCH_SCORES_MULTIBRANCH[basepair_close_inverse.0]
    [basepair_close_inverse.1][basepair_stack_inverse.0][basepair_stack_inverse.1];
  INIT_MULTIBRANCH_BASE
    + terminal_mismatch_score
    + if matches_augu(&basepair_close) {
      HELIX_AUGU_END_PENALTY
    } else {
      0.
    }
}

pub fn get_accessible_score(
  seq: SeqSlice,
  pos_pair_accessible: &(usize, usize),
  uses_sentinel_bases: bool,
) -> Score {
  let seq_len = seq.len();
  let end_5prime = usize::from(uses_sentinel_bases);
  let end_3prime = seq_len - if uses_sentinel_bases { 2 } else { 1 };
  let basepair_accessible = (seq[pos_pair_accessible.0], seq[pos_pair_accessible.1]);
  let score = if pos_pair_accessible.0 > end_5prime && pos_pair_accessible.1 < end_3prime {
    TERMINAL_MISMATCH_SCORES_MULTIBRANCH[basepair_accessible.0][basepair_accessible.1]
      [seq[pos_pair_accessible.0 - 1]][seq[pos_pair_accessible.1 + 1]]
  } else if pos_pair_accessible.0 > end_5prime {
    DANGLING_SCORES_5PRIME[basepair_accessible.0][basepair_accessible.1]
      [seq[pos_pair_accessible.0 - 1]]
  } else if pos_pair_accessible.1 < end_3prime {
    DANGLING_SCORES_3PRIME[basepair_accessible.0][basepair_accessible.1]
      [seq[pos_pair_accessible.1 + 1]]
  } else {
    0.
  };
  score
    + if matches_augu(&basepair_accessible) {
      HELIX_AUGU_END_PENALTY
    } else {
      0.
    }
}

pub fn get_hairpin_score_contra(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  fold_score_sets: &FoldScoreSets,
) -> Score {
  let hairpin_len = pos_pair_close.1 - pos_pair_close.0 - 1;
  fold_score_sets.hairpin_scores_len_cumulative[hairpin_len.min(MAX_LOOP_LEN)]
    + get_junction_score_single(seq, pos_pair_close, fold_score_sets)
}

pub fn get_2loop_score_contra(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
  fold_score_sets: &FoldScoreSets,
) -> Score {
  let basepair_accessible = (seq[pos_pair_accessible.0], seq[pos_pair_accessible.1]);
  let score = if pos_pair_close.0 + 1 == pos_pair_accessible.0
    && pos_pair_close.1 - 1 == pos_pair_accessible.1
  {
    get_stack_score_contra(seq, pos_pair_close, pos_pair_accessible, fold_score_sets)
  } else if pos_pair_close.0 + 1 == pos_pair_accessible.0
    || pos_pair_close.1 - 1 == pos_pair_accessible.1
  {
    get_bulge_score_contra(seq, pos_pair_close, pos_pair_accessible, fold_score_sets)
  } else {
    get_interior_score_contra(seq, pos_pair_close, pos_pair_accessible, fold_score_sets)
  };
  score + fold_score_sets.basepair_scores[basepair_accessible.0][basepair_accessible.1]
}

pub fn get_stack_score_contra(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
  fold_score_sets: &FoldScoreSets,
) -> Score {
  let basepair_close = (seq[pos_pair_close.0], seq[pos_pair_close.1]);
  let basepair_accessible = (seq[pos_pair_accessible.0], seq[pos_pair_accessible.1]);
  fold_score_sets.stack_scores[basepair_close.0][basepair_close.1][basepair_accessible.0]
    [basepair_accessible.1]
}

pub fn get_bulge_score_contra(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
  fold_score_sets: &FoldScoreSets,
) -> Score {
  let bulge_len =
    pos_pair_accessible.0 - pos_pair_close.0 + pos_pair_close.1 - pos_pair_accessible.1 - 2;
  let score = if bulge_len == 1 {
    fold_score_sets.bulge_scores_0x1[if pos_pair_accessible.0 - pos_pair_close.0 - 1 == 1 {
      seq[pos_pair_close.0 + 1]
    } else {
      seq[pos_pair_close.1 - 1]
    }]
  } else {
    0.
  };
  score
    + fold_score_sets.bulge_scores_len_cumulative[bulge_len - 1]
    + get_junction_score_single(seq, pos_pair_close, fold_score_sets)
    + get_junction_score_single(
      seq,
      &(pos_pair_accessible.1, pos_pair_accessible.0),
      fold_score_sets,
    )
}

pub fn get_interior_score_contra(
  seq: SeqSlice,
  pos_pair_close: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
  fold_score_sets: &FoldScoreSets,
) -> Score {
  let loop_len_pair = (
    pos_pair_accessible.0 - pos_pair_close.0 - 1,
    pos_pair_close.1 - pos_pair_accessible.1 - 1,
  );
  let interior_len = loop_len_pair.0 + loop_len_pair.1;
  let score = if loop_len_pair.0 == loop_len_pair.1 {
    let score_1x1 = if interior_len == 2 {
      fold_score_sets.interior_scores_1x1[seq[pos_pair_close.0 + 1]][seq[pos_pair_close.1 - 1]]
    } else {
      0.
    };
    score_1x1 + fold_score_sets.interior_scores_symmetric_cumulative[loop_len_pair.0 - 1]
  } else {
    fold_score_sets.interior_scores_asymmetric_cumulative
      [get_abs_diff(loop_len_pair.0, loop_len_pair.1) - 1]
  };
  let score_explicit =
    if loop_len_pair.0 <= MAX_INTERIOR_EXPLICIT && loop_len_pair.1 <= MAX_INTERIOR_EXPLICIT {
      fold_score_sets.interior_scores_explicit[loop_len_pair.0 - 1][loop_len_pair.1 - 1]
    } else {
      0.
    };
  score
    + score_explicit
    + fold_score_sets.interior_scores_len_cumulative[interior_len - 2]
    + get_junction_score_single(seq, pos_pair_close, fold_score_sets)
    + get_junction_score_single(
      seq,
      &(pos_pair_accessible.1, pos_pair_accessible.0),
      fold_score_sets,
    )
}

pub fn get_junction_score(
  seq: SeqSlice,
  pos_pair: &(usize, usize),
  uses_sentinel_bases: bool,
  fold_score_sets: &FoldScoreSets,
) -> Score {
  let seq_len = seq.len();
  let basepair = (seq[pos_pair.0], seq[pos_pair.1]);
  let end_5prime = usize::from(uses_sentinel_bases);
  let end_3prime = seq_len - if uses_sentinel_bases { 2 } else { 1 };
  get_helix_close_score(&basepair, fold_score_sets)
    + if pos_pair.0 < end_3prime {
      fold_score_sets.dangling_scores_left[basepair.0][basepair.1][seq[pos_pair.0 + 1]]
    } else {
      0.
    }
    + if pos_pair.1 > end_5prime {
      fold_score_sets.dangling_scores_right[basepair.0][basepair.1][seq[pos_pair.1 - 1]]
    } else {
      0.
    }
}

pub fn get_junction_score_single(x: SeqSlice, y: &(usize, usize), z: &FoldScoreSets) -> Score {
  let a = (x[y.0], x[y.1]);
  get_helix_close_score(&a, z) + get_terminal_mismatch_score(&a, &(x[y.0 + 1], x[y.1 - 1]), z)
}

pub fn get_helix_close_score(x: &Basepair, y: &FoldScoreSets) -> Score {
  y.helix_close_scores[x.0][x.1]
}

pub fn get_terminal_mismatch_score(x: &Basepair, y: &Basepair, z: &FoldScoreSets) -> Score {
  z.terminal_mismatch_scores[x.0][x.1][y.0][y.1]
}

pub fn matches_augu(x: &Basepair) -> bool {
  *x == AU || *x == UA || *x == GU || *x == UG
}

pub fn bytes2seq(x: &[Char]) -> Seq {
  let mut y = Seq::new();
  for &x in x {
    let x = match x {
      A_LOWER | A_UPPER => A,
      C_LOWER | C_UPPER => C,
      G_LOWER | G_UPPER => G,
      U_LOWER | U_UPPER => U,
      _ => {
        panic!();
      }
    };
    y.push(x);
  }
  y
}

#[inline]
pub fn logsumexp(sum: &mut Score, x: Score) {
  if !x.is_finite() {
    return;
  }
  *sum = if !sum.is_finite() {
    x
  } else {
    let y = sum.min(x);
    let z = sum.max(x) - y;
    y + if z >= LOGSUMEXP_THRESHOLD_UPPER {
      z
    } else {
      // Equivalent to z.exp().ln_1p()
      ln_exp_1p(z)
    }
  };
}

/*
 * Approximated (x.exp() + 1).ln() from CONTRAfold, eliminating ln() and exp()
 * (assuming 0 <= x <= LOGSUMEXP_THRESHOLD_UPPER)
 */
#[inline]
pub fn ln_exp_1p(x: Score) -> Score {
  if x < 3.379_25 {
    if x < 1.632_015_8 {
      if x < 0.661_536_75 {
        ((-0.0065591595 * x + 0.127_644_27) * x + 0.499_655_46) * x + 0.693_154_2
      } else {
        ((-0.015_515_756 * x + 0.144_677_56) * x + 0.488_293_98) * x + 0.695_809_3
      }
    } else if x < 2.491_258_9 {
      ((-0.012_890_925 * x + 0.130_102_83) * x + 0.515_039_86) * x + 0.679_558_6
    } else {
      ((-0.0072142647 * x + 0.087_754_086) * x + 0.620_870_8) * x + 0.590_967_6
    }
  } else if x < 5.789_071 {
    if x < 4.426_169 {
      ((-0.0031455354 * x + 0.046_722_945) * x + 0.759_253_2) * x + 0.434_879_45
    } else {
      ((-0.0010110698 * x + 0.018_594_341) * x + 0.883_173_05) * x + 0.252_369_55
    }
  } else if x < 7.816_272_7 {
    ((-0.000_196_278 * x + 0.0046084408) * x + 0.963_443_2) * x + 0.098_314_89
  } else {
    ((-0.0000113994 * x + 0.0003734731) * x + 0.995_910_7) * x + 0.0149855051
  }
}

// Approximated x.exp() from CONTRAfold
#[inline]
pub fn expf(x: Score) -> Score {
  if x < -2.491_503_5 {
    if x < -5.862_282_3 {
      if x < -9.91152 {
        0.
      } else {
        ((0.0000803850 * x + 0.002_162_743) * x + 0.019_470_856) * x + 0.058_808_003
      }
    } else if x < -3.839_663 {
      ((0.0013889414 * x + 0.024_467_647) * x + 0.147_129_06) * x + 0.304_275_78
    } else {
      ((0.0072335607 * x + 0.090_600_27) * x + 0.398_311_14) * x + 0.624_595_94
    }
  } else if x < -0.672_505_3 {
    if x < -1.480_537_5 {
      ((0.023_241_036 * x + 0.208_564_6) * x + 0.690_636_8) * x + 0.868_232_25
    } else {
      ((0.057_378_277 * x + 0.358_025_85) * x + 0.912_113_3) * x + 0.979_309_2
    }
  } else if x < 0. {
    ((0.119_917_594 * x + 0.481_566_82) * x + 0.997_599_2) * x + 0.999_950_5
  } else {
    x.exp()
  }
}

pub fn read_align_clustal(clustal_file_path: &Path) -> (Cols, SeqIds) {
  let mut cols = Cols::new();
  let mut seq_ids = SeqIds::new();
  let reader = BufReader::new(File::open(clustal_file_path).unwrap());
  let mut seq_pointer = 0;
  let mut pos_pointer = 0;
  let mut has_read_seq_ids = false;
  for (i, line) in reader.lines().enumerate() {
    let line = line.unwrap();
    if i == 0 || line.is_empty() || line.starts_with(' ') {
      if !cols.is_empty() {
        seq_pointer = 0;
        pos_pointer = cols.len();
        has_read_seq_ids = true;
      }
      continue;
    }
    let mut lines = line.split_whitespace();
    let line = lines.next().unwrap();
    if !has_read_seq_ids {
      seq_ids.push(String::from(line));
    }
    let line = lines.next().unwrap();
    if seq_pointer == 0 {
      for x in line.chars() {
        cols.push(vec![align_char2base(x as Char)]);
      }
      seq_pointer += 1;
    } else {
      for (j, x) in line.chars().enumerate() {
        cols[pos_pointer + j].push(align_char2base(x as Char));
      }
    }
  }
  (cols, seq_ids)
}

pub fn read_align_fasta(fasta_file_path: &Path) -> (Cols, SeqIds) {
  let mut cols = Cols::new();
  let mut seq_ids = SeqIds::new();
  let reader = BufReader::new(File::open(fasta_file_path).unwrap());
  let mut seqs = Vec::<Seq>::new();
  for (i, split) in reader.split(b'>').enumerate() {
    let split = String::from_utf8(split.unwrap()).unwrap();
    if i == 0 {
      continue;
    }
    let splits: Vec<&str> = split.split_whitespace().collect();
    let seq_id = splits[0];
    seq_ids.push(SeqId::from(seq_id));
    let seq = splits[1..].join("");
    let seq = seq.chars().map(|x| align_char2base(x as Char)).collect();
    seqs.push(seq);
  }
  let align_len = seqs[0].len();
  for i in 0..align_len {
    let x = seqs.iter().map(|x| x[i]).collect();
    cols.push(x);
  }
  (cols, seq_ids)
}

pub fn read_align_stockholm(stockholm_file_path: &Path) -> (Cols, SeqIds) {
  let mut cols = Cols::new();
  let mut seq_ids = SeqIds::new();
  let reader = BufReader::new(File::open(stockholm_file_path).unwrap());
  let mut seqs = Vec::<Seq>::new();
  for line in reader.lines() {
    let line = line.unwrap();
    if line.is_empty() || line.starts_with('#') {
      continue;
    } else if line.starts_with("//") {
      break;
    }
    let lines: Vec<&str> = line.split_whitespace().collect();
    let seq_id = lines[0];
    seq_ids.push(SeqId::from(seq_id));
    let seq = lines[1];
    let seq = seq.chars().map(|x| align_char2base(x as Char)).collect();
    seqs.push(seq);
  }
  let align_len = seqs[0].len();
  for i in 0..align_len {
    let x = seqs.iter().map(|x| x[i]).collect();
    cols.push(x);
  }
  (cols, seq_ids)
}

pub fn align_char2base(x: Char) -> Base {
  match x {
    A_LOWER | A_UPPER => A,
    C_LOWER | C_UPPER => C,
    G_LOWER | G_UPPER => G,
    U_LOWER | U_UPPER => U,
    _ => PSEUDO_BASE,
  }
}
