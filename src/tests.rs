extern crate rna_algos;


#[cfg(test)]
mod tests {
  use rna_algos::mccaskill_algo::*;
  
  #[test]
  fn test_mccaskill_algo() {
    let seq = String::from("AUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCC").into_bytes();
    let bpp_matrix = mccaskill_algo(&seq);
    let mut unpaired_base_probs = Probs::new();
    println!("{:?}.", &bpp_matrix);
    for i in 1 .. seq.len() {
      let mut unpaired_base_prob = 1.;
      for j in 1 .. i {
        unpaired_base_prob -= bpp_matrix[j][i];
      }
      for j in i + 1 .. seq.len() {
        unpaired_base_prob -= bpp_matrix[i][j];
      }
      unpaired_base_probs.push(unpaired_base_prob);
    }
    println!("{:?}.", &unpaired_base_probs);
  }
}
