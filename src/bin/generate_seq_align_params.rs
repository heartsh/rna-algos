extern crate rna_algos;

use rna_algos::utils::*;
use std::env;

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "An input file path containing trained CONTRAlign scoring parameters", "STR");
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").unwrap();
  let input_file_path = Path::new(&input_file_path);
  let mut match_score_mat = vec![vec![0.; NUM_OF_BASES]; NUM_OF_BASES];
  let mut insert_scores = vec![0.; NUM_OF_BASES];
  let mut match_2_match_score = 0.;
  let mut match_2_insert_score = 0.;
  let mut insert_extend_score = 0.;
  let mut insert_switch_score = 0.;
  let mut init_match_score = 0.;
  let mut init_insert_score = 0.;
  let reader_2_input_file_path = BufReader::new(File::open(input_file_path).unwrap());
  for line in reader_2_input_file_path.lines() {
    let line = line.unwrap();
    let mut sublines = line.split_whitespace();
    let feature_name = sublines.next().unwrap();
    let feature_weight = sublines.next().unwrap().parse().unwrap();
    if feature_name.eq("match_to_match") {
      match_2_match_score = feature_weight;
    } else if feature_name.eq("match_to_insert") || feature_name.eq("insert_extend") || feature_name.eq("insert_change") || feature_name.eq("insert") {
    } else if feature_name.eq("match_to_insert2") {
      match_2_insert_score = feature_weight;
    } else if feature_name.eq("insert2_extend") {
      insert_extend_score = feature_weight;
    } else if feature_name.eq("insert2_change") {
      insert_switch_score = feature_weight;
    } else if feature_name.eq("match") {
      init_match_score = feature_weight;
    } else if feature_name.eq("insert2") {
      init_insert_score = feature_weight;
    } else if feature_name.starts_with("match_") {
      let prefix = "match_";
      let prefix_len = prefix.len();
      let suffix = feature_name.split_at(prefix_len).1;
      let suffix_chars = suffix.chars().map(|x| x as u8).collect::<Vec<u8>>();
      let char_pair = (convert_char(suffix_chars[0]), convert_char(suffix_chars[1]));
      match_score_mat[char_pair.0][char_pair.1] = feature_weight;
      match_score_mat[char_pair.1][char_pair.0] = feature_weight;
    } else if feature_name.starts_with("insert_") {
      let prefix = "insert_";
      let prefix_len = prefix.len();
      let suffix = feature_name.split_at(prefix_len).1;
      let suffix_chars = suffix.chars().map(|x| x as u8).collect::<Vec<u8>>();
      let base = convert_char(suffix_chars[0]);
      insert_scores[base] = feature_weight;
    } else {
      println!("{}: {}", feature_name, feature_weight);
      assert!(false);
    }
  }
  let output_file_path = Path::new("./src/compiled_seq_align_params.rs");
  let mut writer_2_output_file = BufWriter::new(File::create(&output_file_path).unwrap());
  let mut buf = format!("use utils::*;\n");
  buf += &format!("pub const MATCH_SCORE_MAT: MatchScoreMat = {:?};\n", &match_score_mat);
  buf += &format!("pub const INSERT_SCORES: InsertScores = {:?};\n", &insert_scores);
  buf += &format!("pub const INIT_MATCH_SCORE: Prob = {};\n", &init_match_score);
  buf += &format!("pub const INIT_INSERT_SCORE: Prob = {};\n", init_insert_score);
  buf += &format!("pub const MATCH_2_MATCH_SCORE: Prob = {};\n", match_2_match_score);
  buf += &format!("pub const MATCH_2_INSERT_SCORE: Prob = {};\n", match_2_insert_score);
  buf += &format!("pub const INSERT_EXTEND_SCORE: Prob = {};\n", insert_extend_score);
  buf += &format!("pub const INSERT_SWITCH_SCORE: Prob = {};\n", insert_switch_score);
  let _ = writer_2_output_file.write_all(buf.as_bytes());
}
