extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;

use rna_algos::utils::*;
use rna_algos::mccaskill_algo::*;
use time::precise_time_s;

lazy_static! {
  static ref TEST_SEQ: Seq = {
    convert(&String::from("AUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCC").into_bytes())
  };
  static ref TEST_SEQ_LEN: usize = {TEST_SEQ.len()};
}

#[test]
fn test_mccaskill() {
  let begin = precise_time_s();
  let _ = mccaskill_algo::<u16>(&TEST_SEQ[..], false);
  let elapsed_time = precise_time_s() - begin;
  println!("The elapsed time of McCaskill'a algorithm with Turner's model = {} [s].", elapsed_time);
  let begin = precise_time_s();
  let _ = mccaskill_algo::<u16>(&TEST_SEQ[..], true);
  let elapsed_time = precise_time_s() - begin;
  println!("The elapsed time of McCaskill'a algorithm with CONTRAfold model = {} [s].", elapsed_time);
}
