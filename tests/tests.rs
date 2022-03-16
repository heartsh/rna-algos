extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;

use rna_algos::utils::*;
use rna_algos::mccaskill_algo::*;
use rna_algos::durbin_algo::*;
use time::precise_time_s;

type RealSeqPair = (Seq, Seq);

lazy_static! {
  static ref TEST_SEQ: Seq = {
    convert(&String::from("AUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCC").into_bytes())
  };
  static ref TEST_SEQ_LEN: usize = {TEST_SEQ.len()};
  static ref TEST_SEQ_PAIR: RealSeqPair = {
    let mut seq = convert(&String::from("AUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCC").into_bytes());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    let mut seq_2 = convert(&String::from("CAAGGCGGCUUAAACACUUGGGAUGCAUGCAAGGGCGCUUUGACACAAGGUAUCCAAGCAAGCGGGCUUAAACUCAUGC").into_bytes());
    seq_2.insert(0, PSEUDO_BASE);
    seq_2.push(PSEUDO_BASE);
    (
      seq,
      seq_2,
    )
  };
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

#[test]
fn test_durbin() {
  let begin = precise_time_s();
  let _ = durbin_algo(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]));
  let elapsed_time = precise_time_s() - begin;
  println!("The elapsed time  = {} [s].", elapsed_time);
}
