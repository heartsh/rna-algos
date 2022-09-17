extern crate rna_algos;

pub use rna_algos::durbin_algo::*;
pub use rna_algos::utils::*;

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "An input FASTA file path containing RNA sequences", "STR");
  opts.reqopt("o", "output_file_path", "An output file path", "STR");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
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
  let output_file_path = matches.opt_str("o").unwrap();
  let output_file_path = Path::new(&output_file_path);
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = convert(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut align_feature_score_sets = AlignFeatureCountSets::new(0.);
  align_feature_score_sets.transfer();
  let num_of_fasta_records = fasta_records.len();
  let mut thread_pool = Pool::new(num_of_threads);
  let mut align_prob_mats_with_rna_id_pairs = ProbMatsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      align_prob_mats_with_rna_id_pairs.insert(rna_id_pair, ProbMat::new());
    }
  }
  let ref ref_2_align_feature_score_sets = align_feature_score_sets;
  thread_pool.scoped(|scope| {
    for (rna_id_pair, align_prob_mat) in &mut align_prob_mats_with_rna_id_pairs {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      scope.execute(move || {
        *align_prob_mat = durbin_algo(&seq_pair, ref_2_align_feature_score_sets);
      });
    }
  });
  let mut buf_4_writer_2_align_prob_mat_file = format!("# Format = >{{RNA sequence id 1}},{{RNA sequence id 2}} {{line break}} {{nucleotide 1}}, {{nucleotide 2}}, {{nucletide matching probability}} ...");
  let mut writer_2_align_prob_mat_file = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id_pair, align_prob_mat) in &align_prob_mats_with_rna_id_pairs {
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (i, align_probs) in align_prob_mat.iter().enumerate() {
      for (j, &align_prob) in align_probs.iter().enumerate() {
        if align_prob > 0. {
          buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, align_prob));
        }
      }
    }
    buf_4_writer_2_align_prob_mat_file.push_str(&buf_4_rna_id_pair);
  }
  let _ = writer_2_align_prob_mat_file.write_all(buf_4_writer_2_align_prob_mat_file.as_bytes());
}
