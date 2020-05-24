pub use rna_ss_params::utils::*;
pub use rna_ss_params::hairpin_loop_params::*;
pub use rna_ss_params::helix_params::*;
pub use rna_ss_params::dangling_end_params::*;
pub use rna_ss_params::compiled_free_energy_params::*;
pub use std::cmp::{min, max};
pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use itertools::multizip;

pub type Pos = u16;
pub type PosPair = (Pos, Pos);
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
pub type FastaRecords = Vec<FastaRecord>;
pub type SeqSlice<'a> = &'a[Base];

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

impl FastaRecord {
  pub fn new(fasta_id: FastaId, seq: Seq) -> FastaRecord {
    FastaRecord {
      fasta_id: fasta_id,
      seq: seq,
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

pub fn is_rna_base(base: Base) -> bool {
  match base {
    A => true,
    U => true,
    G => true,
    C => true,
    _ => false,
  }
}

pub fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
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
