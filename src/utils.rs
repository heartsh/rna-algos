use rna_ss_params::utils;
pub use rna_ss_params::utils::*;
pub use rna_ss_params::hairpin_loop_params::*;
use rna_ss_params::terminal_mismatch_params::*;
use rna_ss_params::helix_params::*;
use rna_ss_params::bulge_loop_params::*;
use rna_ss_params::internal_loop_params::*;
pub use bio_seq_algos::utils::*;
pub use std::cmp::{min, max};

pub type Pos = usize;
pub type PosPair = (Pos, Pos);
pub type Num = usize;
type NumPair = (Num, Num);
pub type Energy = Prob;

pub const MAX_IL_LEN: usize = 30;
pub const MAX_SPAN_OF_INDEX_PAIR_CLOSING_IL: usize = MAX_IL_LEN + 2;
pub const INVERSE_TEMPERATURE: Energy = 1. / (GAS_CONST as Energy * TEMPERATURE as Energy); // The unit is [K * mol / (kcal * K)] = [mol / kcal].
pub const MIN_SPAN_OF_INDEX_PAIR_CLOSING_ML: usize = MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL * 2 + 2;
lazy_static! {
  static ref CANONICAL_BPS: HashMap<BasePair, bool> = {
    [(AU, true), (CG, true), (GC, true), (GU, true), (UA, true), (UG, true)].iter().cloned().collect()
  };
}

#[inline]
pub fn is_canonical(bp: &BasePair) -> bool {
  match CANONICAL_BPS.get(bp) {
    Some(_) => true,
    None => false
  }
}

#[inline]
pub fn get_hl_fe(seq: SeqSlice, pp_closing_loop: &PosPair) -> FreeEnergy {
  let hl = &seq[pp_closing_loop.0 .. pp_closing_loop.1 + 1];
  match SPECIAL_HL_DELTA_FES.get(hl) {
    Some(&fe) => fe,
    None => {
      let hl_len = pp_closing_loop.1 - pp_closing_loop.0 - 1;
      if hl_len == MIN_HL_LEN {
        INIT_HL_DELTA_FES[hl_len] + if is_all_c_hl(hl) {HL_OF_3_CS_PENALTY_DELTA_FE} else {0.}
      } else {
        let bp_closing_hl = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
        let tm = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
        let seq_len = seq.len();
        let pair_of_4_bases_preceding_bp_closing_hl = if pp_closing_loop.0 > 1 && pp_closing_loop.1 < seq_len - 2 {
          Some((
            (seq[pp_closing_loop.0 - 2], seq[pp_closing_loop.0 - 1]),
            (seq[pp_closing_loop.1 + 1], seq[pp_closing_loop.1 + 2])
          ))
        } else {None};
        INIT_HL_DELTA_FES[hl_len] 
        + TM_DELTA_FES[&(bp_closing_hl, tm)]
        + if tm == UU || tm == GA {HL_UU_OR_GA_FIRST_MISMATCH_BONUS_DELTA_FE} else if tm == GG {HL_GG_FIRST_MISMATCH_BONUS_DELTA_FE} else {0.}
        + if bp_closing_hl == GU {
          match pair_of_4_bases_preceding_bp_closing_hl {
            Some((bp_1, bp_2)) => if bp_1 == GG || bp_2 == GG {
              HL_SPECIAL_GU_CLOSURE_BONUS_DELTA_FE
            } else {
              0.
            },
            None => 0.
          }
        } else {0.}
        + if is_all_c_hl(hl) {COEFFICENT_4_ALL_C_HL_DELTA_FE * hl_len as FreeEnergy + CONST_4_ALL_C_HL_DELTA_FE} else {0.}
      }
    }
  }
}

#[inline]
pub fn get_2_loop_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  if pp_closing_loop.0 + 1 == accessible_pp.0 && pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_stack_fe(seq, pp_closing_loop, accessible_pp)
  } else if pp_closing_loop.0 + 1 == accessible_pp.0 || pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_bl_fe(seq, pp_closing_loop, accessible_pp)
  } else {
    get_il_fe(seq, pp_closing_loop, accessible_pp)
  }
}

#[inline]
fn get_stack_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  STACK_DELTA_FES[&(bp_closing_loop, accessible_bp)]
}

#[inline]
fn get_bl_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bl_len = accessible_pp.0 - pp_closing_loop.0 + pp_closing_loop.1 - accessible_pp.1 - 2;
  if bl_len == 1 {
    let bulge = if pp_closing_loop.0 + 1 == accessible_pp.0 {seq[accessible_pp.1 + 1]} else {seq[accessible_pp.0 - 1]};
    let bp_adjacent_2_bulge = if pp_closing_loop.0 + 1 == accessible_pp.0 {(seq[accessible_pp.1], seq[pp_closing_loop.1])} else {(seq[pp_closing_loop.0], seq[accessible_pp.0])};
    INIT_BL_DELTA_FES[bl_len]
    + if (pp_closing_loop.0 + 1 == accessible_pp.0 && bulge == C && (bp_adjacent_2_bulge.0 == C || bp_adjacent_2_bulge.1 == C)) || (pp_closing_loop.1 - 1 == accessible_pp.1 && bulge == C && (bp_adjacent_2_bulge.0 == C || bp_adjacent_2_bulge.1 == C)) {BL_SPECIAL_C_BULGE_BONUS_DELTA_FE} else {0.}
    + get_stack_fe(seq, pp_closing_loop, accessible_pp)
  } else {
    if bl_len <= MAX_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_LOOP_DELTA_FE {
      INIT_BL_DELTA_FES[bl_len]
    } else {
      INIT_BL_DELTA_FES[MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_BL_DELTA_FE - 1] + COEFFICENT_4_LOG_EXTRAPOLATION_OF_INIT_BL_DELTA_FE * utils::fast_ln(bl_len as FreeEnergy / (MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_BL_DELTA_FE - 1) as FreeEnergy)
    }
  }
}

#[inline]
fn get_il_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  let pair_of_nums_of_unpaired_bases = (accessible_pp.0 - pp_closing_loop.0 - 1, pp_closing_loop.1 - accessible_pp.1 - 1);
  let il_len = pair_of_nums_of_unpaired_bases.0 + pair_of_nums_of_unpaired_bases.1;
  match il_len {
    2 => {
      let il = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
      ONE_VS_1_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    3 => {
      let il = if pair_of_nums_of_unpaired_bases.0 == 1 {
        ((seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]), seq[pp_closing_loop.1 - 2])
      } else {
        ((seq[pp_closing_loop.1 - 1], seq[pp_closing_loop.0 + 2]), seq[pp_closing_loop.0 + 1])
      };
      ONE_VS_2_IL_DELTA_FES[&if pair_of_nums_of_unpaired_bases.0 == 1 {(bp_closing_loop, il, accessible_bp)} else {(invert_bp(&accessible_bp), il, invert_bp(&bp_closing_loop))}]
    },
    4 => {
      let il = (
        (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
        (seq[pp_closing_loop.0 + 2], seq[pp_closing_loop.1 - 2])
      );
      TWO_VS_2_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    _ => {
      get_init_il_delta_fe(il_len)
      + IL_ASYMMETRY_PENALTY_DELTA_FE * abs(pair_of_nums_of_unpaired_bases.0, pair_of_nums_of_unpaired_bases.1) as FreeEnergy
      + get_il_tm_bonus_delta_fe(seq, pp_closing_loop, accessible_pp, &pair_of_nums_of_unpaired_bases)
      + if bp_closing_loop == AU || bp_closing_loop == GU {IL_AU_OR_GU_CLOSURE_PENALTY_DELTA_FE} else {0.}
      + if accessible_bp == AU || accessible_bp == GU {IL_AU_OR_GU_CLOSURE_PENALTY_DELTA_FE} else {0.}
    },
  }
}

#[inline]
fn get_init_il_delta_fe(il_len: usize) -> FreeEnergy {
  if il_len <= MAX_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_LOOP_DELTA_FE {
    INIT_IL_DELTA_FES[il_len]
  } else {
    INIT_IL_DELTA_FES[MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_IL_DELTA_FE - 1] + COEFFICENT_4_LOG_EXTRAPOLATION_OF_INIT_IL_DELTA_FE * utils::fast_ln(il_len as FreeEnergy / (MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_IL_DELTA_FE - 1) as FreeEnergy)
  }
}

#[inline]
fn invert_bp(bp: &BasePair) -> BasePair {(bp.1, bp.0)}

#[inline]
fn abs(x: usize, y: usize) -> usize {max(x, y) - min(x, y)}

#[inline]
fn get_il_tm_bonus_delta_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair, pair_of_nums_of_unpaired_bases: &NumPair) -> FreeEnergy {
  let tm_pair = (
    (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
    (seq[accessible_pp.0 - 1], seq[accessible_pp.1 + 1])
  );
  match *pair_of_nums_of_unpaired_bases {
    (1, _) => 0.,
    (_, 1) => 0.,
    (2, 3) => {
      match TWO_VS_3_IL_TM_BONUS_DELTA_FES.get(&tm_pair.0) {
        Some(&fe_1) => {
          match TWO_VS_3_IL_TM_BONUS_DELTA_FES.get(&tm_pair.1) {
            Some(&fe_2) => fe_1 + fe_2,
            None => fe_1 
          }
        },
        None => 0.
      }
    },
    _ => {
      match OTHER_IL_TM_BONUS_DELTA_FES.get(&tm_pair.0) {
        Some(&fe_1) => {
          match OTHER_IL_TM_BONUS_DELTA_FES.get(&tm_pair.1) {
            Some(&fe_2) => fe_1 + fe_2,
            None => fe_1
          }
        },
        None => 0.
      }
    }
  }
}

#[inline]
fn is_all_c_hl(hl: SeqSlice) -> bool {
  for base in hl {if *base != C {return false;}}
  return true;
}
