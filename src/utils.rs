pub use rna_ss_params::utils::*;
pub use rna_ss_params::compiled_free_energy_params_turner::*;
pub use rna_ss_params::compiled_free_energy_params_contra::*;
pub use std::cmp::{min, max};
pub use std::str::from_utf8_unchecked;
pub use itertools::multizip;
pub use num::{Unsigned, PrimInt, One, Zero, FromPrimitive, ToPrimitive, Bounded, range, range_inclusive, Integer};
pub use std::hash::Hash;
pub use std::fmt::Display;
pub use scoped_threadpool::Pool;
pub use std::env;
pub use std::path::Path;
pub use bio::io::fasta::Reader;
pub use std::fs::File;
pub use std::fs::create_dir;
pub use hashbrown::HashMap;
pub use std::f32::NEG_INFINITY;
pub use std::io::prelude::*;
pub use std::io::{BufReader, BufWriter};

pub type PosPair<T> = (T, T);
pub type PosQuadruple<T> = (T, T, T, T);
pub type Num = usize;
type NumPair = (Num, Num);
type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
#[derive(Clone)]
pub struct FastaRecord {
  pub fasta_id: FastaId,
  pub seq: Seq,
}
#[derive(Debug)]
pub struct SeqAlign<T> {
  pub cols: Cols,
  pub pos_map_sets: PosMapSets<T>,
}
pub type PosMaps<T> = Vec<T>;
pub type PosMapSets<T> = Vec<PosMaps<T>>;
pub type FastaRecords = Vec<FastaRecord>;
pub type SeqSlice<'a> = &'a[Base];
pub type NumOfThreads = u32;
pub type SparseProbMat<T> = HashMap<PosPair<T>, Prob>;
pub type BaScoreMat = HashMap<BasePair, FreeEnergy>;
pub type BpaScoreMat = HashMap<(BasePair, BasePair), FreeEnergy>;
pub type SeqId = String;
pub type SeqIds = Vec<SeqId>;
pub type Col = Vec<Base>;
pub type Cols = Vec<Col>;
pub type PartFunc4dMat<T> = HashMap<PosQuadruple<T>, PartFunc>;
pub type SparsePartFuncMat<T> = HashMap<PosPair<T>, PartFunc>;
pub type FreeEnergies = Vec<FreeEnergy>;
pub type FreeEnergyMat = Vec<FreeEnergies>;
pub type SparseFreeEnergyMat<T> = HashMap<PosPair<T>, FreeEnergy>;
pub type PosPairs<T> = Vec<PosPair<T>>;
pub type Mea = Prob;
pub type MeaSsChar = u8;
pub type MeaSsStr = Vec<MeaSsChar>;
pub type Char = u8;
pub type Seq = Vec<Base>;
pub type SeqPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type Prob = f32;
pub type LogProb = Prob;
pub type PartFunc = Prob;
pub type Probs = Vec<Prob>;
pub type ProbMat = Vec<Probs>;
pub type PartFuncs = Vec<PartFunc>;
pub type PartFuncMat = Vec<PartFuncs>;
pub type Pos = usize;
pub type Base = usize;
pub type MatchScoreMat = [[Prob; NUM_OF_BASES]; NUM_OF_BASES];
pub type InsertScores = [Prob; NUM_OF_BASES];
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type ProbMatsWithRnaIdPairs = HashMap<RnaIdPair, ProbMat>;

pub const MAX_SPAN_OF_INDEX_PAIR_CLOSING_IL: usize = MAX_2_LOOP_LEN + 2;
pub const MIN_SPAN_OF_INDEX_PAIR_CLOSING_ML: usize = MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL * 2 + 2;
pub const SMALL_A: u8 = 'a' as u8;
pub const BIG_A: u8 = 'A' as u8;
pub const SMALL_C: u8 = 'c' as u8;
pub const BIG_C: u8 = 'C' as u8;
pub const SMALL_G: u8 = 'g' as u8;
pub const BIG_G: u8 = 'G' as u8;
pub const SMALL_U: u8 = 'u' as u8;
pub const BIG_U: u8 = 'U' as u8;
pub const LOGSUMEXP_THRES_UPPER: FreeEnergy = 11.8624794162;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const U: Base = 3;
pub const PSEUDO_BASE: Base = U + 1 as Base;
pub const UNPAIRING_BASE: MeaSsChar = '.' as MeaSsChar;
pub const BASE_PAIRING_LEFT_BASE: MeaSsChar = '(' as MeaSsChar;
pub const BASE_PAIRING_RIGHT_BASE: MeaSsChar = ')' as MeaSsChar;
pub const NUM_OF_BASES: usize = 4;

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

impl<T> SeqAlign<T> {
  pub fn new() -> SeqAlign<T> {
    SeqAlign {
      cols: Cols::new(),
      pos_map_sets: PosMapSets::<T>::new(),
    }
  }
}

pub fn is_canonical(bp: &BasePair) -> bool {
  match *bp {
    AU | CG | GC | GU | UA | UG => true,
    _ => false,
  }
}

pub fn get_hl_fe(seq: SeqSlice, pp_closing_loop: &(usize, usize)) -> FreeEnergy {
  let hl = &seq[pp_closing_loop.0 .. pp_closing_loop.1 + 1];
  let special_hl_fe = get_special_hl_fe(hl);
  if special_hl_fe > NEG_INFINITY {
    special_hl_fe
  } else {
    let hl_len = pp_closing_loop.1 - pp_closing_loop.0 - 1;
    let bp_closing_hl = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
    let hl_fe = if hl_len == MIN_HL_LEN {
      INIT_HL_DELTA_FES[hl_len]
    } else {
      let tm = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
      let init_hl_delta_fe = if hl_len <= MAX_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_LOOP_DELTA_FE {
        INIT_HL_DELTA_FES[hl_len]
      } else {
        INIT_HL_DELTA_FES[MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_HL_DELTA_FE - 1] + COEFFICIENT_4_LOG_EXTRAPOLATION_OF_INIT_HL_DELTA_FE * (hl_len as FreeEnergy / (MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_HL_DELTA_FE - 1) as FreeEnergy).ln()
      };
      init_hl_delta_fe + HL_TM_DELTA_FES[bp_closing_hl.0][bp_closing_hl.1][tm.0][tm.1]
    };
    hl_fe + if is_au_or_gu(&bp_closing_hl) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
  }
}

pub fn get_special_hl_fe(seq: SeqSlice) -> FreeEnergy {
  for special_hl_delta_fe in SPECIAL_HL_DELTA_FES.iter() {
    if special_hl_delta_fe.0 == seq {
      return special_hl_delta_fe.1;
    }
  }
  NEG_INFINITY
}

pub fn get_2_loop_fe(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  if pp_closing_loop.0 + 1 == accessible_pp.0 && pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_stack_fe(seq, pp_closing_loop, accessible_pp)
  } else if pp_closing_loop.0 + 1 == accessible_pp.0 || pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_bl_fe(seq, pp_closing_loop, accessible_pp)
  } else {
    get_il_fe(seq, pp_closing_loop, accessible_pp)
  }
}

fn get_stack_fe(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  STACK_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][accessible_bp.0][accessible_bp.1]
}

fn get_bl_fe(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let bl_len = accessible_pp.0 - pp_closing_loop.0 + pp_closing_loop.1 - accessible_pp.1 - 2;
  if bl_len == 1 {
    INIT_BL_DELTA_FES[bl_len]
    + get_stack_fe(seq, pp_closing_loop, accessible_pp)
  } else {
    let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
    let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
    INIT_BL_DELTA_FES[bl_len] + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
    + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
  }
}

fn get_il_fe(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  let pair_of_nums_of_unpaired_bases = (accessible_pp.0 - pp_closing_loop.0 - 1, pp_closing_loop.1 - accessible_pp.1 - 1);
  let il_len = pair_of_nums_of_unpaired_bases.0 + pair_of_nums_of_unpaired_bases.1;
  match pair_of_nums_of_unpaired_bases {
    (1, 1) => {
      let il = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
      ONE_VS_1_IL_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][il.0][il.1][accessible_bp.0][accessible_bp.1]
    },
    (1, 2) => {
      let il = ((seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]), seq[pp_closing_loop.1 - 2]);
      ONE_VS_2_IL_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][(il.0).0][(il.0).1][il.1][accessible_bp.0][accessible_bp.1]
    },
    (2, 1) => {
      let il = ((seq[pp_closing_loop.1 - 1], seq[pp_closing_loop.0 + 2]), seq[pp_closing_loop.0 + 1]);
      let invert_accessible_bp = invert_bp(&accessible_bp);
      let invert_bp_closing_loop = invert_bp(&bp_closing_loop);
      ONE_VS_2_IL_DELTA_FES[invert_accessible_bp.0][invert_accessible_bp.1][(il.0).0][(il.0).1][il.1][invert_bp_closing_loop.0][invert_bp_closing_loop.1]
    },
    (2, 2) => {
      let il = (
        (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
        (seq[pp_closing_loop.0 + 2], seq[pp_closing_loop.1 - 2])
      );
      TWO_VS_2_IL_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][(il.0).0][(il.0).1][(il.1).0][(il.1).1][accessible_bp.0][accessible_bp.1]
    },
    _ => {
      INIT_IL_DELTA_FES[il_len]
      + (COEFFICIENT_4_NINIO * get_abs_diff(pair_of_nums_of_unpaired_bases.0, pair_of_nums_of_unpaired_bases.1) as FreeEnergy).max(MAX_NINIO)
      + get_il_tm_delta_fe(seq, pp_closing_loop, accessible_pp, &pair_of_nums_of_unpaired_bases)
      + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
      + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
    },
  }
}

pub fn invert_bp(bp: &BasePair) -> BasePair {(bp.1, bp.0)}

pub fn get_abs_diff(x: usize, y: usize) -> usize {max(x, y) - min(x, y)}

fn get_il_tm_delta_fe(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize), pair_of_nums_of_unpaired_bases: &NumPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.1], seq[accessible_pp.0]);
  let tm_pair = (
    (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
    (seq[accessible_pp.1 + 1], seq[accessible_pp.0 - 1]),
  );
  match *pair_of_nums_of_unpaired_bases {
    (1, _) => {
      ONE_VS_MANY_IL_TM_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][(tm_pair.0).0][(tm_pair.0).1] + ONE_VS_MANY_IL_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][(tm_pair.1).0][(tm_pair.1).1]
    },
    (_, 1) => {
      ONE_VS_MANY_IL_TM_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][(tm_pair.0).0][(tm_pair.0).1] + ONE_VS_MANY_IL_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][(tm_pair.1).0][(tm_pair.1).1]
    },
    (2, 3) => {
      TWO_VS_3_IL_TM_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][(tm_pair.0).0][(tm_pair.0).1] + TWO_VS_3_IL_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][(tm_pair.1).0][(tm_pair.1).1]
    },
    (3, 2) => {
      TWO_VS_3_IL_TM_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][(tm_pair.0).0][(tm_pair.0).1] + TWO_VS_3_IL_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][(tm_pair.1).0][(tm_pair.1).1]
    },
    _ => {
      IL_TM_DELTA_FES[bp_closing_loop.0][bp_closing_loop.1][(tm_pair.0).0][(tm_pair.0).1] + IL_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][(tm_pair.1).0][(tm_pair.1).1]
    }
  }
}

pub fn get_ml_closing_basepairing_fe(seq: SeqSlice, pp_closing_loop: &(usize, usize)) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let invert_bp_closing_loop = invert_bp(&bp_closing_loop);
  let invert_stacking_bp = invert_bp(&(seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]));
  let ml_tm_delta_fe = ML_TM_DELTA_FES[invert_bp_closing_loop.0][invert_bp_closing_loop.1][invert_stacking_bp.0][invert_stacking_bp.1];
  CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe + if is_au_or_gu(&bp_closing_loop) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
}

pub fn get_ml_or_el_accessible_basepairing_fe(seq: SeqSlice, pp_accessible: &(usize, usize), uses_sentinel_nucs: bool) -> FreeEnergy {
  let seq_len = seq.len();
  let five_prime_end = if uses_sentinel_nucs {1} else {0};
  let three_prime_end = seq_len - if uses_sentinel_nucs {2} else {1};
  let accessible_bp = (seq[pp_accessible.0], seq[pp_accessible.1]);
  let fe = if pp_accessible.0 > five_prime_end && pp_accessible.1 < three_prime_end {
    ML_TM_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[pp_accessible.0 - 1]][seq[pp_accessible.1 + 1]]
  } else if pp_accessible.0 > five_prime_end {
    FIVE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[pp_accessible.0 - 1]]
  } else if pp_accessible.1 < three_prime_end {
    THREE_PRIME_DE_DELTA_FES[accessible_bp.0][accessible_bp.1][seq[pp_accessible.1 + 1]]
  } else {
    0.
  } + if is_au_or_gu(&accessible_bp) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
  fe
}

pub fn get_hl_fe_contra(seq: SeqSlice, pp_closing_loop: &(usize, usize)) -> FreeEnergy {
  let hl_len = pp_closing_loop.1 - pp_closing_loop.0 - 1;
  CONTRA_HL_LENGTH_FES[hl_len.min(CONTRA_MAX_LOOP_LEN)]
    + get_contra_junction_fe_single(seq, pp_closing_loop)
}

pub fn get_2_loop_fe_contra(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  let fe = if pp_closing_loop.0 + 1 == accessible_pp.0 && pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_stack_fe_contra(seq, pp_closing_loop, accessible_pp)
  } else if pp_closing_loop.0 + 1 == accessible_pp.0 || pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_bl_fe_contra(seq, pp_closing_loop, accessible_pp)
  } else {
    get_il_fe_contra(seq, pp_closing_loop, accessible_pp)
  };
  fe + CONTRA_BASE_PAIR_FES[accessible_bp.0][accessible_bp.1]
}

pub fn get_stack_fe_contra(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  CONTRA_STACK_FES[bp_closing_loop.0][bp_closing_loop.1][accessible_bp.0][accessible_bp.1]
}

pub fn get_bl_fe_contra(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let bl_len = accessible_pp.0 - pp_closing_loop.0 + pp_closing_loop.1 - accessible_pp.1 - 2;
  let fe = if bl_len == 1 {
    CONTRA_BL_0X1_FES[if accessible_pp.0 - pp_closing_loop.0 - 1 == 1 {
      seq[pp_closing_loop.0 + 1]
    } else {
      seq[pp_closing_loop.1 - 1]
    }]
  } else {0.};
  fe + CONTRA_BL_LENGTH_FES[bl_len - 1]
    + get_contra_junction_fe_single(seq, pp_closing_loop)
    + get_contra_junction_fe_single(seq, &(accessible_pp.1, accessible_pp.0))
}

pub fn get_il_fe_contra(seq: SeqSlice, pp_closing_loop: &(usize, usize), accessible_pp: &(usize, usize)) -> FreeEnergy {
  let pair_of_nums_of_unpaired_bases = (accessible_pp.0 - pp_closing_loop.0 - 1, pp_closing_loop.1 - accessible_pp.1 - 1);
  let il_len = pair_of_nums_of_unpaired_bases.0 + pair_of_nums_of_unpaired_bases.1;
  let fe = if pair_of_nums_of_unpaired_bases.0 == pair_of_nums_of_unpaired_bases.1 {
    let fe_3 = if il_len == 2 {
      CONTRA_IL_1X1_FES[seq[pp_closing_loop.0 + 1]][seq[pp_closing_loop.1 - 1]]
    } else {0.};
    fe_3 + CONTRA_IL_SYMM_LENGTH_FES[pair_of_nums_of_unpaired_bases.0 - 1]
  } else {
    CONTRA_IL_ASYMM_LENGTH_FES[get_abs_diff(pair_of_nums_of_unpaired_bases.0, pair_of_nums_of_unpaired_bases.1) - 1]
  };
  let fe_2 = if pair_of_nums_of_unpaired_bases.0 <= 4 && pair_of_nums_of_unpaired_bases.1 <= 4 {
    CONTRA_IL_EXPLICIT_FES[pair_of_nums_of_unpaired_bases.0 - 1][pair_of_nums_of_unpaired_bases.1 - 1]
  } else {0.};
  fe + fe_2 + CONTRA_IL_LENGTH_FES[il_len - 2]
    + get_contra_junction_fe_single(seq, pp_closing_loop)
    + get_contra_junction_fe_single(seq, &(accessible_pp.1, accessible_pp.0))
}

pub fn get_contra_junction_fe_multi(seq: SeqSlice, pp: &(usize, usize), seq_len: usize, uses_sentinel_nucs: bool) -> FreeEnergy {
  let bp = (seq[pp.0], seq[pp.1]);
  let five_prime_end = if uses_sentinel_nucs {1} else {0};
  let three_prime_end = seq_len - if uses_sentinel_nucs {2} else {1};
  get_contra_helix_closing_fe(&bp) + if pp.0 < three_prime_end && pp.1 > five_prime_end {
    CONTRA_LEFT_DANGLE_FES[bp.0][bp.1][seq[pp.0 + 1]] + CONTRA_RIGHT_DANGLE_FES[bp.0][bp.1][seq[pp.1 - 1]]
  } else if pp.0 < three_prime_end {
    CONTRA_LEFT_DANGLE_FES[bp.0][bp.1][seq[pp.0 + 1]]
  } else if pp.1 > five_prime_end {
    CONTRA_RIGHT_DANGLE_FES[bp.0][bp.1][seq[pp.1 - 1]]
  } else {
    0.
  }
}

pub fn get_contra_junction_fe_single(seq: SeqSlice, pp: &(usize, usize)) -> FreeEnergy {
  let bp = (seq[pp.0], seq[pp.1]);
  get_contra_helix_closing_fe(&bp) + get_contra_terminal_mismatch_fe(&bp, &(seq[pp.0 + 1], seq[pp.1 - 1]))
}

pub fn get_contra_helix_closing_fe(bp: &BasePair) -> FreeEnergy {
  CONTRA_HELIX_CLOSING_FES[bp.0][bp.1]
}

pub fn get_contra_terminal_mismatch_fe(bp: &BasePair, mismatch_bp: &BasePair) -> FreeEnergy {
  CONTRA_TERMINAL_MISMATCH_FES[bp.0][bp.1][mismatch_bp.0][mismatch_bp.1]
}

pub fn is_rna_base(base: Base) -> bool {
  match base {
    A => true,
    U => true,
    G => true,
    C => true,
    _ => false,
  }
}

pub fn is_au_or_gu(bp: &BasePair) -> bool {
  *bp == AU || *bp == UA || *bp == GU || *bp == UG
}

pub fn convert<'a>(seq: &'a [u8]) -> Seq {
  let mut new_seq = Seq::new();
  for &c in seq {
    let new_base = match c {
      SMALL_A | BIG_A => A,
      SMALL_C | BIG_C => C,
      SMALL_G | BIG_G => G,
      SMALL_U | BIG_U => U,
      _ => {assert!(false); U},
    };
    new_seq.push(new_base);
  }
  new_seq
}

#[inline]
pub fn logsumexp(sum: &mut FreeEnergy, new_term: FreeEnergy) {
  if !new_term.is_finite() {
    return;
  }
  *sum = if !sum.is_finite() {
    new_term
  } else {
    let max = sum.max(new_term);
    let min = sum.min(new_term);
    let diff = max - min;
    min + if diff >= LOGSUMEXP_THRES_UPPER {
      diff
    } else {
      // diff.exp().ln_1p()
      ln_exp_1p(diff)
    }
  };
}

// Approximated (x.exp() + 1).ln() from CONTRAfold, eliminating ln() and exp() (assuming 0 <= x <= LOGSUMEXP_THRES_UPPER)
#[inline]
pub fn ln_exp_1p(x: FreeEnergy) -> FreeEnergy {
  if x < 3.3792499610 {
    if x < 1.6320158198 {
      if x < 0.6615367791 {
        ((-0.0065591595 * x + 0.1276442762) * x + 0.4996554598) * x + 0.6931542306
      } else {
        ((-0.0155157557 * x + 0.1446775699) * x + 0.4882939746) * x + 0.6958092989
      }
    } else if x < 2.4912588184 {
      ((-0.0128909247 * x + 0.1301028251) * x + 0.5150398748) * x + 0.6795585882
    } else {
      ((-0.0072142647 * x + 0.0877540853) * x + 0.6208708362) * x + 0.5909675829
    }
  } else if x < 5.7890710412 {
    if x < 4.4261691294 {
      ((-0.0031455354 * x + 0.0467229449) * x + 0.7592532310) * x + 0.4348794399
    } else {
      ((-0.0010110698 * x + 0.0185943421) * x + 0.8831730747) * x + 0.2523695427
    }
  } else if x < 7.8162726752 {
    ((-0.0001962780 * x + 0.0046084408) * x + 0.9634431978) * x + 0.0983148903
  } else {
    ((-0.0000113994 * x + 0.0003734731) * x + 0.9959107193) * x + 0.0149855051
  }
}

// Approximated x.exp() from CONTRAfold
#[inline]
pub fn expf(x: FreeEnergy) -> FreeEnergy {
  if x < -2.4915033807 {
    if x < -5.8622823336 {
      if x < -9.91152 {
        0.
      } else {
        ((0.0000803850 * x + 0.0021627428) * x + 0.0194708555) * x + 0.0588080014
      }
    } else if x < -3.8396630909 {
      ((0.0013889414 * x + 0.0244676474) * x + 0.1471290604) * x + 0.3042757740
    } else {
      ((0.0072335607 * x + 0.0906002677) * x + 0.3983111356) * x + 0.6245959221
    }
  } else if x < -0.6725053211 {
    if x < -1.4805375919 {
      ((0.0232410351 * x + 0.2085645908) * x + 0.6906367911) * x + 0.8682322329
    } else {
      ((0.0573782771 * x + 0.3580258429) * x + 0.9121133217) * x + 0.9793091728
    }
  } else if x < 0. {
    ((0.1199175927 * x + 0.4815668234) * x + 0.9975991939) * x + 0.9999505077
  } else {
    x.exp()
  }
}

pub fn read_sa_from_clustal_file(clustal_file_path: &Path) -> (Cols, SeqIds) {
  let mut cols = Cols::new();
  let mut seq_ids = SeqIds::new();
  let reader_2_clustal_file = BufReader::new(File::open(clustal_file_path).unwrap());
  let mut seq_pointer = 0;
  let mut pos_pointer = 0;
  let mut are_seq_ids_read = false;
  for (i, string) in reader_2_clustal_file.lines().enumerate() {
    let string = string.unwrap();
    if i == 0 || string.len() == 0 || string.starts_with(" ") {
      if cols.len() > 0 {
        seq_pointer = 0;
        pos_pointer = cols.len();
        are_seq_ids_read = true;
      }
      continue;
    }
    let mut substrings = string.split_whitespace();
    let substring = substrings.next().unwrap();
    if !are_seq_ids_read {
      seq_ids.push(String::from(substring));
    }
    let substring = substrings.next().unwrap();
    if seq_pointer == 0 {
      for sa_char in substring.chars() {
        cols.push(vec![convert_sa_char(sa_char as u8)]);
      }
      seq_pointer += 1;
    } else {
      for (j, sa_char) in substring.chars().enumerate() {
        cols[pos_pointer + j].push(convert_sa_char(sa_char as u8));
      }
    }
  }
  (cols, seq_ids)
}

pub fn read_sa_from_fasta_file(fasta_file_path: &Path) -> (Cols, SeqIds) {
  let mut cols = Cols::new();
  let mut seq_ids = SeqIds::new();
  let reader_2_fasta_file = BufReader::new(File::open(fasta_file_path).unwrap());
  let mut seqs = Vec::<Seq>::new();
  for (i, string) in reader_2_fasta_file.split(b'>').enumerate() {
    let string = String::from_utf8(string.unwrap()).unwrap();
    if i == 0 {continue;}
    let substrings: Vec<&str> = string.split_whitespace().collect();
    let seq_id = substrings[0];
    seq_ids.push(SeqId::from(seq_id));
    let seq = substrings[1 ..].join("");
    let seq = seq.chars().map(|x| convert_sa_char(x as u8)).collect();
    seqs.push(seq);
  }
  let align_len = seqs[0].len();
  for i in 0 .. align_len {
    let col = seqs.iter().map(|x| x[i]).collect();
    cols.push(col);
  }
  (cols, seq_ids)
}

pub fn read_sa_from_stockholm_file(stockholm_file_path: &Path) -> (Cols, SeqIds) {
  let mut cols = Cols::new();
  let mut seq_ids = SeqIds::new();
  let reader_2_stockholm_file = BufReader::new(File::open(stockholm_file_path).unwrap());
  let mut seqs = Vec::<Seq>::new();
  for string in reader_2_stockholm_file.lines() {
    let string = string.unwrap();
    if string.len() == 0 || string.starts_with("#") {
      continue;
    } else if string.starts_with("//") {
      break;
    }
    let substrings: Vec<&str> = string.split_whitespace().collect();
    let seq_id = substrings[0];
    seq_ids.push(SeqId::from(seq_id));
    let seq = substrings[1];
    let seq = seq.chars().map(|x| convert_sa_char(x as u8)).collect();
    seqs.push(seq);
  }
  let align_len = seqs[0].len();
  for i in 0 .. align_len {
    let col = seqs.iter().map(|x| x[i]).collect();
    cols.push(col);
  }
  (cols, seq_ids)
}

pub fn convert_sa_char(c: u8) -> Base {
  match c {
    SMALL_A | BIG_A => A,
    SMALL_C | BIG_C => C,
    SMALL_G | BIG_G => G,
    SMALL_U | BIG_U => U,
    _ => {PSEUDO_BASE},
  }
}
