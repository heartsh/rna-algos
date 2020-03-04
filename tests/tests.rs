extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;

use rna_algos::utils::*;
use rna_algos::mccaskill_algo::*;
use time::precise_time_s;

lazy_static! {
  static ref TEST_SEQ: Seq = {
    convert(&String::from("AUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCC").into_bytes())
  };
  static ref TEST_SEQ_LEN: usize = {TEST_SEQ.len()};
}

#[test]
fn test_bpp_and_upp_mats() {
  let begin = precise_time_s();
  let _ = get_bpp_and_unpair_prob_mats(&TEST_SEQ[..]);
  let elapsed_time = precise_time_s() - begin;
  println!("The elapsed time = {} [s].", elapsed_time);
}
