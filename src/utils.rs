pub use rna_ss_params::utils::*;
pub use rna_ss_params::hairpin_loop_params::*;
pub use rna_ss_params::helix_params::*;
use rna_ss_params::bulge_loop_params::*;
pub use rna_ss_params::dangling_end_params::*;
use rna_ss_params::internal_loop_params::*;
pub use std::cmp::{min, max};
pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use itertools::multizip;

pub type Pos = usize;
pub type PosPair = (Pos, Pos);
pub type Num = usize;
type NumPair = (Num, Num);
pub type HlTmDeltaFes = StackDeltaFes;
pub type IlTmDeltaFes = HlTmDeltaFes;
pub type MlTmDeltaFes = HlTmDeltaFes;
type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
#[derive(Clone)]
pub struct FastaRecord {
  pub fasta_id: FastaId,
  pub seq: Seq,
}
pub type FastaRecords = Vec<FastaRecord>;

pub const MAX_SPAN_OF_INDEX_PAIR_CLOSING_IL: usize = MAX_2_LOOP_LEN + 2;
pub const MIN_SPAN_OF_INDEX_PAIR_CLOSING_ML: usize = MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL * 2 + 2;
pub const HELIX_AU_OR_GU_END_PENALTY_DELTA_FE: FreeEnergy = - INVERSE_TEMPERATURE * 0.5;
pub const MAX_NINIO: FreeEnergy = - INVERSE_TEMPERATURE * 3.;
pub const COEFFICIENT_4_NINIO: FreeEnergy = - INVERSE_TEMPERATURE * 0.6;
lazy_static! {
  pub static ref EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE: FreeEnergy = HELIX_AU_OR_GU_END_PENALTY_DELTA_FE.exp();
  pub static ref EXP_MAX_NINIO: FreeEnergy = MAX_NINIO.exp();
  pub static ref EXP_COEFFICIENT_4_NINIO: FreeEnergy = COEFFICIENT_4_NINIO.exp();
  static ref CANONICAL_BPS: HashMap<BasePair, bool> = {
    [(AU, true), (CG, true), (GC, true), (GU, true), (UA, true), (UG, true)].iter().cloned().collect()
  };
  pub static ref HL_TM_DELTA_FES: HlTmDeltaFes = {
    [
      // For the base pair "AU" against which another base pair is stacked.
      ((AU, AA), -0.3),  ((AU, AC), -0.5),  ((AU, AG), -0.3),  ((AU, AU), -0.5),
      ((AU, CA), -0.1),  ((AU, CC), -0.2),  ((AU, CG), -0.1),  ((AU, CU), -0.2),
      ((AU, GA), -1.2),  ((AU, GC), -0.5), ((AU, GG), -1.1),  ((AU, GU), -0.5),
      ((AU, UA), -0.1),  ((AU, UC), -0.3),  ((AU, UG), -0.1), ((AU, UU), -1.2),
      // For the base pair "CG" against which another base pair is stacked.
      ((CG, AA), -1.5), ((CG, AC), -1.5), ((CG, AG), -1.4), ((CG, AU), -1.5),
      ((CG, CA), -1.0), ((CG, CC), -1.1), ((CG, CG), -1.0), ((CG, CU), -0.8),
      ((CG, GA), -2.3), ((CG, GC), -1.5), ((CG, GG), -2.4), ((CG, GU), -1.5),
      ((CG, UA), -1.0), ((CG, UC), -1.4), ((CG, UG), -1.0), ((CG, UU), -2.1),
      // For the base pair "GC" against which another base pair is stacked.
      ((GC, AA), -1.1), ((GC, AC), -1.5), ((GC, AG), -1.3), ((GC, AU), -1.5),
      ((GC, CA), -1.1),  ((GC, CC), -0.7), ((GC, CG), -1.1),  ((GC, CU), -0.5),
      ((GC, GA), -2.5), ((GC, GC), -1.5), ((GC, GG), -2.2), ((GC, GU), -1.5),
      ((GC, UA), -1.1), ((GC, UC), -1.0), ((GC, UG), -1.1), ((GC, UU), -1.6),
      // For the base pair "GU" against which another base pair is stacked.
      ((GU, AA), 0.2),  ((GU, AC), -0.5),  ((GU, AG), -0.3),  ((GU, AU), -0.5),
      ((GU, CA), -0.1),  ((GU, CC), -0.2),  ((GU, CG), -0.1),  ((GU, CU), -0.2),
      ((GU, GA), -1.0),  ((GU, GC), -0.5), ((GU, GG), -1.1),  ((GU, GU), -0.5),
      ((GU, UA), -0.1),  ((GU, UC), -0.3),  ((GU, UG), -0.1), ((GU, UU), -1.0),
      // For the base pair "UG" against which another base pair is stacked.
      ((UG, AA), -0.5),  ((UG, AC), -0.3),  ((UG, AG), -0.6),  ((UG, AU), -0.3),
      ((UG, CA), -0.2),  ((UG, CC), -0.1),  ((UG, CG), -0.2),    ((UG, CU), 0.0),
      ((UG, GA), -0.9),  ((UG, GC), -0.3), ((UG, GG), -1.1),  ((UG, GU), -0.3),
      ((UG, UA), -0.2),  ((UG, UC), -0.1),  ((UG, UG), -0.2),  ((UG, UU), -0.9),
      // For the base pair "UA" against which another base pair is stacked.
      ((UA, AA), -0.5),  ((UA, AC), -0.3),  ((UA, AG), -0.5),  ((UA, AU), -0.3),
      ((UA, CA), -0.2),  ((UA, CC), -0.1),  ((UA, CG), -0.2),    ((UA, CU), 0.0),
      ((UA, GA), -1.5),  ((UA, GC), -0.3), ((UA, GG), -1.5),  ((UA, GU), -0.3),
      ((UA, UA), -0.2),  ((UA, UC), -0.1),  ((UA, UG), -0.2),  ((UA, UU), -0.9),
    ].iter().map(|&(x, y)| {(x, scale(y))}).collect()
  };
  pub static ref EXP_HL_TM_DELTA_FES: HlTmDeltaFes = {HL_TM_DELTA_FES.iter().map(|(x, &y)| {(*x, y.exp())}).collect()};
  pub static ref IL_TM_DELTA_FES: IlTmDeltaFes = {
    [
      // For the base pair "AU" against which another base pair is stacked.
      ((AU, AA), 0.7),  ((AU, AC), 0.7),  ((AU, AG), -0.1),  ((AU, AU), 0.7),
      ((AU, CA), 0.7),  ((AU, CC), 0.7),  ((AU, CG), 0.7),  ((AU, CU), 0.7),
      ((AU, GA), -0.3),  ((AU, GC), 0.7), ((AU, GG), -0.3),  ((AU, GU), 0.7),
      ((AU, UA), 0.7),  ((AU, UC), 0.7),  ((AU, UG), 0.7), ((AU, UU), 0.1),
      // For the base pair "CG" against which another base pair is stacked.
      ((CG, AA), 0.0), ((CG, AC), 0.0), ((CG, AG), -0.8), ((CG, AU), 0.0),
      ((CG, CA), 0.0), ((CG, CC), 0.0), ((CG, CG), 0.0), ((CG, CU), 0.0),
      ((CG, GA), -1.0), ((CG, GC), 0.0), ((CG, GG), -1.0), ((CG, GU), 0.0),
      ((CG, UA), 0.0), ((CG, UC), 0.0), ((CG, UG), 0.0), ((CG, UU), -0.6),
      // For the base pair "GC" against which another base pair is stacked.
      ((GC, AA), 0.0), ((GC, AC), 0.0), ((GC, AG), -0.8), ((GC, AU), 0.0),
      ((GC, CA), 0.0),  ((GC, CC), 0.0), ((GC, CG), 0.0),  ((GC, CU), 0.0),
      ((GC, GA), -1.0), ((GC, GC), 0.0), ((GC, GG), -1.0), ((GC, GU), 0.0),
      ((GC, UA), 0.0), ((GC, UC), 0.0), ((GC, UG), 0.0), ((GC, UU), -0.6),
      // For the base pair "GU" against which another base pair is stacked.
      ((GU, AA), 0.7),  ((GU, AC), 0.7),  ((GU, AG), -0.1),  ((GU, AU), 0.7),
      ((GU, CA), 0.7),  ((GU, CC), 0.7),  ((GU, CG), 0.7),  ((GU, CU), 0.7),
      ((GU, GA), -0.3),  ((GU, GC), 0.7), ((GU, GG), -0.3),  ((GU, GU), 0.7),
      ((GU, UA), 0.7),  ((GU, UC), 0.7),  ((GU, UG), 0.7), ((GU, UU), 0.1),
      // For the base pair "UG" against which another base pair is stacked.
      ((UG, AA), 0.7),  ((UG, AC), 0.7),  ((UG, AG), -0.1),  ((UG, AU), 0.7),
      ((UG, CA), 0.7),  ((UG, CC), 0.7),  ((UG, CG), 0.7),    ((UG, CU), 0.7),
      ((UG, GA), -0.3),  ((UG, GC), 0.7), ((UG, GG), -0.3),  ((UG, GU), 0.7),
      ((UG, UA), 0.7),  ((UG, UC), 0.7),  ((UG, UG), 0.7),  ((UG, UU), 0.1),
      // For the base pair "UA" against which another base pair is stacked.
      ((UA, AA), 0.7),  ((UA, AC), 0.7),  ((UA, AG), -0.1),  ((UA, AU), 0.7),
      ((UA, CA), 0.7),  ((UA, CC), 0.7),  ((UA, CG), 0.7),    ((UA, CU), 0.7),
      ((UA, GA), -0.3),  ((UA, GC), 0.7), ((UA, GG), -0.3),  ((UA, GU), 0.7),
      ((UA, UA), 0.7),  ((UA, UC), 0.7),  ((UA, UG), 0.7),  ((UA, UU), 0.1),
    ].iter().map(|&(x, y)| {(x, scale(y))}).collect()
  };
  pub static ref EXP_IL_TM_DELTA_FES: IlTmDeltaFes = {IL_TM_DELTA_FES.iter().map(|(x, &y)| {(*x, y.exp())}).collect()};
  pub static ref ONE_VS_MANY_IL_TM_DELTA_FES: IlTmDeltaFes = {
    [
      // For the base pair "AU" against which another base pair is stacked.
      ((AU, AA), 0.7),  ((AU, AC), 0.7),  ((AU, AG), 0.7),  ((AU, AU), 0.7),
      ((AU, CA), 0.7),  ((AU, CC), 0.7),  ((AU, CG), 0.7),  ((AU, CU), 0.7),
      ((AU, GA), 0.7),  ((AU, GC), 0.7), ((AU, GG), 0.7),  ((AU, GU), 0.7),
      ((AU, UA), 0.7),  ((AU, UC), 0.7),  ((AU, UG), 0.7), ((AU, UU), 0.7),
      // For the base pair "CG" against which another base pair is stacked.
      ((CG, AA), 0.0), ((CG, AC), 0.0), ((CG, AG), 0.0), ((CG, AU), 0.0),
      ((CG, CA), 0.0), ((CG, CC), 0.0), ((CG, CG), 0.0), ((CG, CU), 0.0),
      ((CG, GA), 0.0), ((CG, GC), 0.0), ((CG, GG), 0.0), ((CG, GU), 0.0),
      ((CG, UA), 0.0), ((CG, UC), 0.0), ((CG, UG), 0.0), ((CG, UU), 0.0),
      // For the base pair "GC" against which another base pair is stacked.
      ((GC, AA), 0.0), ((GC, AC), 0.0), ((GC, AG), 0.0), ((GC, AU), 0.0),
      ((GC, CA), 0.0),  ((GC, CC), 0.0), ((GC, CG), 0.0),  ((GC, CU), 0.0),
      ((GC, GA), 0.0), ((GC, GC), 0.0), ((GC, GG), 0.0), ((GC, GU), 0.0),
      ((GC, UA), 0.0), ((GC, UC), 0.0), ((GC, UG), 0.0), ((GC, UU), 0.0),
      // For the base pair "GU" against which another base pair is stacked.
      ((GU, AA), 0.7),  ((GU, AC), 0.7),  ((GU, AG), 0.7),  ((GU, AU), 0.7),
      ((GU, CA), 0.7),  ((GU, CC), 0.7),  ((GU, CG), 0.7),  ((GU, CU), 0.7),
      ((GU, GA), 0.7),  ((GU, GC), 0.7), ((GU, GG), 0.7),  ((GU, GU), 0.7),
      ((GU, UA), 0.7),  ((GU, UC), 0.7),  ((GU, UG), 0.7), ((GU, UU), 0.7),
      // For the base pair "UG" against which another base pair is stacked.
      ((UG, AA), 0.7),  ((UG, AC), 0.7),  ((UG, AG), 0.7),  ((UG, AU), 0.7),
      ((UG, CA), 0.7),  ((UG, CC), 0.7),  ((UG, CG), 0.7),    ((UG, CU), 0.7),
      ((UG, GA), 0.7),  ((UG, GC), 0.7), ((UG, GG), 0.7),  ((UG, GU), 0.7),
      ((UG, UA), 0.7),  ((UG, UC), 0.7),  ((UG, UG), 0.7),  ((UG, UU), 0.7),
      // For the base pair "UA" against which another base pair is stacked.
      ((UA, AA), 0.7),  ((UA, AC), 0.7),  ((UA, AG), 0.7),  ((UA, AU), 0.7),
      ((UA, CA), 0.7),  ((UA, CC), 0.7),  ((UA, CG), 0.7),    ((UA, CU), 0.7),
      ((UA, GA), 0.7),  ((UA, GC), 0.7), ((UA, GG), 0.7),  ((UA, GU), 0.7),
      ((UA, UA), 0.7),  ((UA, UC), 0.7),  ((UA, UG), 0.7),  ((UA, UU), 0.7),
    ].iter().map(|&(x, y)| {(x, scale(y))}).collect()
  };
  pub static ref EXP_ONE_VS_MANY_IL_TM_DELTA_FES: IlTmDeltaFes = {ONE_VS_MANY_IL_TM_DELTA_FES.iter().map(|(x, &y)| {(*x, y.exp())}).collect()};
  pub static ref TWO_VS_3_IL_TM_DELTA_FES: IlTmDeltaFes = {
    [
      // For the base pair "AU" against which another base pair is stacked.
      ((AU, AA), 0.7),  ((AU, AC), 0.7),  ((AU, AG), 0.7),  ((AU, AU), 0.7),
      ((AU, CA), 0.7),  ((AU, CC), 0.7),  ((AU, CG), 0.7),  ((AU, CU), 0.7),
      ((AU, GA), 0.4),  ((AU, GC), 0.7), ((AU, GG), 0.0),  ((AU, GU), 0.7),
      ((AU, UA), 0.7),  ((AU, UC), 0.7),  ((AU, UG), 0.7), ((AU, UU), 0.4),
      // For the base pair "CG" against which another base pair is stacked.
      ((CG, AA), 0.0), ((CG, AC), 0.0), ((CG, AG), -0.5), ((CG, AU), 0.0),
      ((CG, CA), 0.0), ((CG, CC), 0.0), ((CG, CG), 0.0), ((CG, CU), 0.0),
      ((CG, GA), -1.1), ((CG, GC), 0.0), ((CG, GG), 0.7), ((CG, GU), 0.0),
      ((CG, UA), 0.0), ((CG, UC), 0.0), ((CG, UG), 0.0), ((CG, UU), -0.3),
      // For the base pair "GC" against which another base pair is stacked.
      ((GC, AA), 0.0), ((GC, AC), 0.0), ((GC, AG), 0.0), ((GC, AU), 0.0),
      ((GC, CA), 0.0),  ((GC, CC), 0.0), ((GC, CG), 0.0),  ((GC, CU), 0.0),
      ((GC, GA), -1.2), ((GC, GC), 0.0), ((GC, GG), -0.7), ((GC, GU), 0.0),
      ((GC, UA), 0.0), ((GC, UC), 0.0), ((GC, UG), 0.0), ((GC, UU), -0.3),
      // For the base pair "GU" against which another base pair is stacked.
      ((GU, AA), 0.7),  ((GU, AC), 0.7),  ((GU, AG), 0.7),  ((GU, AU), 0.7),
      ((GU, CA), 0.7),  ((GU, CC), 0.7),  ((GU, CG), 0.7),  ((GU, CU), 0.7),
      ((GU, GA), -0.4),  ((GU, GC), 0.7), ((GU, GG), 0.0),  ((GU, GU), 0.7),
      ((GU, UA), 0.7),  ((GU, UC), 0.7),  ((GU, UG), 0.7), ((GU, UU), 0.4),
      // For the base pair "UG" against which another base pair is stacked.
      ((UG, AA), 0.7),  ((UG, AC), 0.7),  ((UG, AG), 0.7),  ((UG, AU), 0.7),
      ((UG, CA), 0.7),  ((UG, CC), 0.7),  ((UG, CG), 0.7),    ((UG, CU), 0.7),
      ((UG, GA), -0.4),  ((UG, GC), 0.7), ((UG, GG), 0.0),  ((UG, GU), 0.7),
      ((UG, UA), 0.7),  ((UG, UC), 0.7),  ((UG, UG), 0.7),  ((UG, UU), 0.4),
      // For the base pair "UA" against which another base pair is stacked.
      ((UA, AA), 0.7),  ((UA, AC), 0.7),  ((UA, AG), 0.7),  ((UA, AU), 0.7),
      ((UA, CA), 0.7),  ((UA, CC), 0.7),  ((UA, CG), 0.7),    ((UA, CU), 0.7),
      ((UA, GA), -0.4),  ((UA, GC), 0.7), ((UA, GG), 0.0),  ((UA, GU), 0.7),
      ((UA, UA), 0.7),  ((UA, UC), 0.7),  ((UA, UG), 0.7),  ((UA, UU), 0.4),
    ].iter().map(|&(x, y)| {(x, scale(y))}).collect()
  };
  pub static ref EXP_TWO_VS_3_IL_TM_DELTA_FES: IlTmDeltaFes = {TWO_VS_3_IL_TM_DELTA_FES.iter().map(|(x, &y)| {(*x, y.exp())}).collect()};
  pub static ref ML_TM_DELTA_FES: MlTmDeltaFes = {
    [
      // For the base pair "AU" against which another base pair is stacked.
      ((AU, AA), -0.8),  ((AU, AC), -1.0),  ((AU, AG), -0.8),  ((AU, AU), -1.0),
      ((AU, CA), -0.6),  ((AU, CC), -0.7),  ((AU, CG), -0.6),  ((AU, CU), -0.7),
      ((AU, GA), -0.8),  ((AU, GC), -1.0), ((AU, GG), -0.8),  ((AU, GU), -1.0),
      ((AU, UA), -0.6),  ((AU, UC), -0.8),  ((AU, UG), -0.6), ((AU, UU), -0.8),
      // For the base pair "CG" against which another base pair is stacked.
      ((CG, AA), -1.5), ((CG, AC), -1.5), ((CG, AG), -1.4), ((CG, AU), -1.5),
      ((CG, CA), -1.0), ((CG, CC), -1.1), ((CG, CG), -1.0), ((CG, CU), -0.8),
      ((CG, GA), -1.4), ((CG, GC), -1.5), ((CG, GG), -1.6), ((CG, GU), -1.5),
      ((CG, UA), -1.0), ((CG, UC), -1.4), ((CG, UG), -1.0), ((CG, UU), -1.2),
      // For the base pair "GC" against which another base pair is stacked.
      ((GC, AA), -1.1), ((GC, AC), -1.5), ((GC, AG), -1.3), ((GC, AU), -1.5),
      ((GC, CA), -1.1),  ((GC, CC), -0.7), ((GC, CG), -1.1),  ((GC, CU), -0.5),
      ((GC, GA), -1.6), ((GC, GC), -1.5), ((GC, GG), -1.4), ((GC, GU), -1.5),
      ((GC, UA), -1.1), ((GC, UC), -1.0), ((GC, UG), -1.1), ((GC, UU), -0.7),
      // For the base pair "GU" against which another base pair is stacked.
      ((GU, AA), -0.3),  ((GU, AC), -1.0),  ((GU, AG), -0.8),  ((GU, AU), -1.0),
      ((GU, CA), -0.6),  ((GU, CC), -0.7),  ((GU, CG), -0.6),  ((GU, CU), -0.7),
      ((GU, GA), -0.6),  ((GU, GC), -1.0), ((GU, GG), -0.8),  ((GU, GU), -1.0),
      ((GU, UA), -0.6),  ((GU, UC), -0.8),  ((GU, UG), -0.6), ((GU, UU), -0.6),
      // For the base pair "UG" against which another base pair is stacked.
      ((UG, AA), -1.0),  ((UG, AC), -0.8),  ((UG, AG), -1.1),  ((UG, AU), -0.8),
      ((UG, CA), -0.7),  ((UG, CC), -0.6),  ((UG, CG), -0.7),    ((UG, CU), -0.5),
      ((UG, GA), -0.5),  ((UG, GC), -0.8), ((UG, GG), -0.8),  ((UG, GU), -0.8),
      ((UG, UA), -0.7),  ((UG, UC), -0.6),  ((UG, UG), -0.7),  ((UG, UU), -0.5),
      // For the base pair "UA" against which another base pair is stacked.
      ((UA, AA), -1.0),  ((UA, AC), -0.8),  ((UA, AG), -1.1),  ((UA, AU), -0.8),
      ((UA, CA), -0.7),  ((UA, CC), -0.6),  ((UA, CG), -0.7),    ((UA, CU), -0.5),
      ((UA, GA), -1.1),  ((UA, GC), -0.8), ((UA, GG), -1.2),  ((UA, GU), -0.8),
      ((UA, UA), -0.7),  ((UA, UC), -0.6),  ((UA, UG), -0.7),  ((UA, UU), -0.5),
    ].iter().map(|&(x, y)| {(x, scale(y))}).collect()
  };
  pub static ref EXP_ML_TM_DELTA_FES: MlTmDeltaFes = {ML_TM_DELTA_FES.iter().map(|(x, &y)| {(*x, y.exp())}).collect()};
}

impl FastaRecord {
  pub fn new(fasta_id: FastaId, seq: Seq) -> FastaRecord {
    FastaRecord {
      fasta_id: fasta_id,
      seq: seq,
    }
  }
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
        init_hl_delta_fe + HL_TM_DELTA_FES[&(bp_closing_hl, tm)]
      };
      hl_fe + if bp_closing_hl == AU || bp_closing_hl == UA || bp_closing_hl == GU || bp_closing_hl == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
    }
  }
}

#[inline]
pub fn get_exp_hl_fe(seq: SeqSlice, pp_closing_loop: &PosPair) -> FreeEnergy {
  let hl = &seq[pp_closing_loop.0 .. pp_closing_loop.1 + 1];
  match EXP_SPECIAL_HL_DELTA_FES.get(hl) {
    Some(&fe) => fe,
    None => {
      let hl_len = pp_closing_loop.1 - pp_closing_loop.0 - 1;
      let bp_closing_hl = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
      let exp_hl_fe = if hl_len == MIN_HL_LEN {
        EXP_INIT_HL_DELTA_FES[hl_len]
      } else {
        let tm = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
        let exp_init_hl_delta_fe = if hl_len <= MAX_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_LOOP_DELTA_FE {
          EXP_INIT_HL_DELTA_FES[hl_len]
        } else {
          EXP_INIT_HL_DELTA_FES[MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_HL_DELTA_FE - 1] * (hl_len as FreeEnergy / (MIN_LOOP_LEN_4_LOG_EXTRAPOLATION_OF_INIT_HL_DELTA_FE - 1) as FreeEnergy).powf(COEFFICIENT_4_LOG_EXTRAPOLATION_OF_INIT_HL_DELTA_FE)
        };
        exp_init_hl_delta_fe * EXP_HL_TM_DELTA_FES[&(bp_closing_hl, tm)]
      };
      exp_hl_fe * if bp_closing_hl == AU || bp_closing_hl == UA || bp_closing_hl == GU || bp_closing_hl == UG {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.}
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
pub fn get_exp_2_loop_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  if pp_closing_loop.0 + 1 == accessible_pp.0 && pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_exp_stack_fe(seq, pp_closing_loop, accessible_pp)
  } else if pp_closing_loop.0 + 1 == accessible_pp.0 || pp_closing_loop.1 - 1 == accessible_pp.1 {
    get_exp_bl_fe(seq, pp_closing_loop, accessible_pp)
  } else {
    get_exp_il_fe(seq, pp_closing_loop, accessible_pp)
  }
}

#[inline]
fn get_stack_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  STACK_DELTA_FES[&(bp_closing_loop, accessible_bp)]
}

#[inline]
fn get_exp_stack_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  EXP_STACK_DELTA_FES[&(bp_closing_loop, accessible_bp)]
}

#[inline]
fn get_bl_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bl_len = accessible_pp.0 - pp_closing_loop.0 + pp_closing_loop.1 - accessible_pp.1 - 2;
  if bl_len == 1 {
    INIT_BL_DELTA_FES[bl_len]
    + get_stack_fe(seq, pp_closing_loop, accessible_pp)
  } else {
    let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
    let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
    INIT_BL_DELTA_FES[bl_len] + if bp_closing_loop == AU || bp_closing_loop == UA || bp_closing_loop == GU || bp_closing_loop == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
    + if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
  }
}

#[inline]
fn get_exp_bl_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bl_len = accessible_pp.0 - pp_closing_loop.0 + pp_closing_loop.1 - accessible_pp.1 - 2;
  if bl_len == 1 {
    EXP_INIT_BL_DELTA_FES[bl_len]
    * get_exp_stack_fe(seq, pp_closing_loop, accessible_pp)
  } else {
    let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
    let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
    EXP_INIT_BL_DELTA_FES[bl_len] * if bp_closing_loop == AU || bp_closing_loop == UA || bp_closing_loop == GU || bp_closing_loop == UG {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.}
    * if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.}
  }
}

#[inline]
fn get_il_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  let pair_of_nums_of_unpaired_bases = (accessible_pp.0 - pp_closing_loop.0 - 1, pp_closing_loop.1 - accessible_pp.1 - 1);
  let il_len = pair_of_nums_of_unpaired_bases.0 + pair_of_nums_of_unpaired_bases.1;
  match pair_of_nums_of_unpaired_bases {
    (1,1) => {
      let il = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
      ONE_VS_1_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    (1, 2) => {
      let il = ((seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]), seq[pp_closing_loop.1 - 2]);
      ONE_VS_2_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    (2, 1) => {
      let il = ((seq[pp_closing_loop.1 - 1], seq[pp_closing_loop.0 + 2]), seq[pp_closing_loop.0 + 1]);
      ONE_VS_2_IL_DELTA_FES[&(invert_bp(&accessible_bp), il, invert_bp(&bp_closing_loop))]
    },
    (2, 2) => {
      let il = (
        (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
        (seq[pp_closing_loop.0 + 2], seq[pp_closing_loop.1 - 2])
      );
      TWO_VS_2_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    _ => {
      INIT_IL_DELTA_FES[il_len]
      + (COEFFICIENT_4_NINIO * get_abs_diff(pair_of_nums_of_unpaired_bases.0, pair_of_nums_of_unpaired_bases.1) as FreeEnergy).min(MAX_NINIO)
      + get_il_tm_delta_fe(seq, pp_closing_loop, accessible_pp, &pair_of_nums_of_unpaired_bases)
      + if bp_closing_loop == AU || bp_closing_loop == UA || bp_closing_loop == GU || bp_closing_loop == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
      + if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
    },
  }
}

#[inline]
fn get_exp_il_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.0], seq[accessible_pp.1]);
  let pair_of_nums_of_unpaired_bases = (accessible_pp.0 - pp_closing_loop.0 - 1, pp_closing_loop.1 - accessible_pp.1 - 1);
  let il_len = pair_of_nums_of_unpaired_bases.0 + pair_of_nums_of_unpaired_bases.1;
  match pair_of_nums_of_unpaired_bases {
    (1,1) => {
      let il = (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]);
      EXP_ONE_VS_1_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    (1, 2) => {
      let il = ((seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]), seq[pp_closing_loop.1 - 2]);
      EXP_ONE_VS_2_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    (2, 1) => {
      let il = ((seq[pp_closing_loop.1 - 1], seq[pp_closing_loop.0 + 2]), seq[pp_closing_loop.0 + 1]);
      EXP_ONE_VS_2_IL_DELTA_FES[&(invert_bp(&accessible_bp), il, invert_bp(&bp_closing_loop))]
    },
    (2, 2) => {
      let il = (
        (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
        (seq[pp_closing_loop.0 + 2], seq[pp_closing_loop.1 - 2])
      );
      EXP_TWO_VS_2_IL_DELTA_FES[&(bp_closing_loop, il, accessible_bp)]
    },
    _ => {
      EXP_INIT_IL_DELTA_FES[il_len]
      * EXP_COEFFICIENT_4_NINIO.powi(get_abs_diff(pair_of_nums_of_unpaired_bases.0, pair_of_nums_of_unpaired_bases.1) as i32).min(*EXP_MAX_NINIO)
      * get_exp_il_tm_delta_fe(seq, pp_closing_loop, accessible_pp, &pair_of_nums_of_unpaired_bases)
      * if bp_closing_loop == AU || bp_closing_loop == UA || bp_closing_loop == GU || bp_closing_loop == UG {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.}
      * if accessible_bp == AU || accessible_bp == UA || accessible_bp == GU || accessible_bp == UG {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.}
    },
  }
}

#[inline]
pub fn invert_bp(bp: &BasePair) -> BasePair {(bp.1, bp.0)}

#[inline]
pub fn get_abs_diff(x: usize, y: usize) -> usize {max(x, y) - min(x, y)}

#[inline]
fn get_il_tm_delta_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair, pair_of_nums_of_unpaired_bases: &NumPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.1], seq[accessible_pp.0]);
  let tm_pair = (
    (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
    (seq[accessible_pp.1 + 1], seq[accessible_pp.0 - 1]),
  );
  match *pair_of_nums_of_unpaired_bases {
    (1, _) => {
      ONE_VS_MANY_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] + ONE_VS_MANY_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    (_, 1) => {
      ONE_VS_MANY_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] + ONE_VS_MANY_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    (2, 3) => {
      TWO_VS_3_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] + TWO_VS_3_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    (3, 2) => {
      TWO_VS_3_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] + TWO_VS_3_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    _ => {
      IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] + IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    }
  }
}

#[inline]
fn get_exp_il_tm_delta_fe(seq: SeqSlice, pp_closing_loop: &PosPair, accessible_pp: &PosPair, pair_of_nums_of_unpaired_bases: &NumPair) -> FreeEnergy {
  let bp_closing_loop = (seq[pp_closing_loop.0], seq[pp_closing_loop.1]);
  let accessible_bp = (seq[accessible_pp.1], seq[accessible_pp.0]);
  let tm_pair = (
    (seq[pp_closing_loop.0 + 1], seq[pp_closing_loop.1 - 1]),
    (seq[accessible_pp.1 + 1], seq[accessible_pp.0 - 1]),
  );
  match *pair_of_nums_of_unpaired_bases {
    (1, _) => {
      EXP_ONE_VS_MANY_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] * EXP_ONE_VS_MANY_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    (_, 1) => {
      EXP_ONE_VS_MANY_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] * EXP_ONE_VS_MANY_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    (2, 3) => {
      EXP_TWO_VS_3_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] * EXP_TWO_VS_3_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    (3, 2) => {
      EXP_TWO_VS_3_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] * EXP_TWO_VS_3_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    },
    _ => {
      EXP_IL_TM_DELTA_FES[&(bp_closing_loop, tm_pair.0)] * EXP_IL_TM_DELTA_FES[&(accessible_bp, tm_pair.1)]
    }
  }
}

#[inline]
pub fn is_rna_base(base: Base) -> bool {
  match base {
    A => true,
    U => true,
    G => true,
    C => true,
    _ => false,
  }
}

#[inline]
pub fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
}

#[inline]
pub fn is_au_or_gu(bp: &BasePair) -> bool {
  *bp == AU || *bp == UA || *bp == GU || *bp == UG
}
