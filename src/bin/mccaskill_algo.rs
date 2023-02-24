extern crate rna_algos;

pub use rna_algos::mccaskill_algo::*;
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
  opts.optflag(
    "c",
    "uses_contra_model",
    "Use the CONTRAfold model instead of Turner's model to score RNA secondary structures",
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
  let uses_contra_model = matches.opt_present("c");
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let seq = bytes2seq(fasta_record.seq());
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let num_fasta_records = fasta_records.len();
  let mut thread_pool = Pool::new(num_threads);
  let mut fold_score_sets = FoldScoreSets::new(0.);
  fold_score_sets.transfer();
  let ref_fold_score_sets = &fold_score_sets;
  let mut basepair_prob_mats = vec![String::new(); num_fasta_records];
  let allows_short_hairpins = false;
  thread_pool.scoped(|scope| {
    for (basepair_probs, fasta_record) in multizip((basepair_prob_mats.iter_mut(), fasta_records.iter())) {
      scope.execute(move || {
        let seq_len = fasta_record.seq.len();
        if seq_len <= u8::MAX as usize {
          *basepair_probs = probs2str(&mccaskill_algo::<u8>(
            &fasta_record.seq[..],
            uses_contra_model,
            allows_short_hairpins,
            ref_fold_score_sets,
          ).0);
        } else {
          *basepair_probs = probs2str(&mccaskill_algo::<u16>(
            &fasta_record.seq[..],
            uses_contra_model,
            allows_short_hairpins,
            ref_fold_score_sets,
          ).0);
        }
      });
    }
  });
  let mut buf = "# Format = >{RNA sequence id} {line break} {basepairing left nucleotide}, {basepairing right nucleotide}, {basepairing probability} ...".to_string();
  let mut writer = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id, basepair_probs) in basepair_prob_mats.iter().enumerate() {
    let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
    buf_rna_id.push_str(basepair_probs);
    buf.push_str(&buf_rna_id);
  }
  let _ = writer.write_all(buf.as_bytes());
}

fn probs2str<T>(probs: &SparseProbMat<T>) -> String
where
  T: HashIndex,
{
  let mut probs_str = String::new();
  for (&(i, j), &x) in probs.iter() {
    probs_str.push_str(&format!("{},{},{} ", &i.to_string(), &j.to_string(), x));
  }
  probs_str
}
