extern crate rna_algos;
#[macro_use]
extern crate lazy_static;

use rna_algos::utils::*;
use rna_algos::mccaskill_algo::*;

lazy_static! {
  static ref TEST_SEQ: Seq = {
    String::from("AUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCC").into_bytes()
  };
  static ref TEST_SEQ_LEN: usize = {TEST_SEQ.len()};
}

#[test]
fn test_bpp_and_upp_mats() {
  let (bpp_mat, upp_mat, _) = get_bpp_and_unpair_prob_mats(&TEST_SEQ[..]);
  println!("The base-pairing matrix for the sequence \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ[..]), &bpp_mat);
  for &bpp in bpp_mat.values() {
    assert!(0. <= bpp && bpp <= 1.);
  }
  println!("The unpairing probabilities for the sequence \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ[..]), &upp_mat);
  for upp in upp_mat {assert!((0. <= upp && upp <= 1.));}
}
