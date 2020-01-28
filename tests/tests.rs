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

/* #[test]
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
} */

#[test]
fn test_bpp_and_upp_mats() {
  let (bpp_mat, upp_mat, _) = get_bpp_and_unpair_prob_mats(&TEST_SEQ[..]);
  println!("The base-pairing matrix for the sequence \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ[..]), &bpp_mat);
  for &bpp in bpp_mat.values() {
    assert!(0. <= bpp && bpp <= 1.);
    // for &bpp in bpps {assert!((0. <= bpp && bpp <= 1.));}
  }
  println!("The not-base-pairing probabilities for the sequence \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ[..]), &upp_mat);
  for upp in upp_mat {assert!((0. <= upp && upp <= 1.));}
}
