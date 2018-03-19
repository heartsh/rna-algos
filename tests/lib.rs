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
fn test_log_ss_ppf_mat() {
  let log_ss_ppf_mat = get_log_ss_ppf_mats(&TEST_SEQ[..], *TEST_SEQ_LEN).log_ss_ppf_mat;
  for i in 0 .. *TEST_SEQ_LEN {
    for j in i + 1 .. *TEST_SEQ_LEN {
      let log_ss_ppf_1 = log_ss_ppf_mat[i][j];
      for k in i + 1 .. j - 1 {
        for l in k + 1 .. j {
          let log_ss_ppf_2 = log_ss_ppf_mat[k][l];
          assert!(log_ss_ppf_1 >= log_ss_ppf_2 || !log_ss_ppf_1.is_finite() || !log_ss_ppf_2.is_finite());
        }
      }
    }
  }
}

#[test]
fn test_bpp_mat_and_ubps() {
  let (bpp_mat, ubps) = get_bpp_mat_and_unpaired_base_probs(&TEST_SEQ[..]);
  println!("The base pairing mat for the seq. \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ[..]), &bpp_mat);
  for bpps in &bpp_mat {
    for &bpp in bpps {assert!((0. <= bpp && bpp <= 1.));}
  }
  println!("The unpaired base probs. for the seq. \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ[..]), &ubps);
  for ubp in ubps {assert!((0. <= ubp && ubp <= 1.));}
}
