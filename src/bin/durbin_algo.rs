extern crate rna_algos;

pub use rna_algos::durbin_algo::*;
pub use rna_algos::utils::*;

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt(
    "i",
    "input_file_path",
    "An input FASTA file path containing RNA sequences",
    "STR",
  );
  opts.reqopt("o", "output_file_path", "An output file path", "STR");
  opts.optopt(
    "t",
    "num_threads",
    "The number of threads in multithreading (Use all the threads of this computer by default)",
    "UINT",
  );
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1..]) {
    Ok(opt) => opt,
    Err(failure) => {
      print_program_usage(&program_name, &opts);
      panic!("{}", failure.to_string())
    }
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").unwrap();
  let input_file_path = Path::new(&input_file_path);
  let output_file_path = matches.opt_str("o").unwrap();
  let output_file_path = Path::new(&output_file_path);
  let num_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = bytes2seq(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut align_scores = AlignScores::new(0.);
  align_scores.transfer();
  let num_fasta_records = fasta_records.len();
  let mut thread_pool = Pool::new(num_threads);
  let mut match_probs_hashed_ids = ProbMatsHashedIds::default();
  for rna_id in 0..num_fasta_records {
    for rna_id2 in rna_id + 1..num_fasta_records {
      let rna_id_pair = (rna_id, rna_id2);
      match_probs_hashed_ids.insert(rna_id_pair, ProbMat::new());
    }
  }
  let ref_align_scores = &align_scores;
  thread_pool.scoped(|scope| {
    for (rna_id_pair, match_prob_mat) in &mut match_probs_hashed_ids {
      let seq_pair = (
        &fasta_records[rna_id_pair.0].seq[..],
        &fasta_records[rna_id_pair.1].seq[..],
      );
      scope.execute(move || {
        *match_prob_mat = durbin_algo(&seq_pair, ref_align_scores);
      });
    }
  });
  let mut buf = "# Format = >{RNA sequence id 1},{RNA sequence id 2} {line break} {nucleotide 1}, {nucleotide 2}, {nucletide matching probability} ...".to_string();
  let mut writer = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id_pair, match_probs) in &match_probs_hashed_ids {
    let mut buf_rna_id = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (i, x) in match_probs.iter().enumerate() {
      for (j, &x) in x.iter().enumerate() {
        if x > 0. {
          buf_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, x));
        }
      }
    }
    buf.push_str(&buf_rna_id);
  }
  let _ = writer.write_all(buf.as_bytes());
}
