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
    "num_of_threads",
    "The number of threads in multithreading (Use all the threads of this computer by default)",
    "UINT",
  );
  opts.optflag(
    "c",
    "use_contra_model",
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
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let use_contra_model = matches.opt_present("c");
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let seq = convert(fasta_record.seq());
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let num_of_fasta_records = fasta_records.len();
  let mut thread_pool = Pool::new(num_of_threads);
  let mut struct_feature_score_sets = StructFeatureCountSets::new(0.);
  struct_feature_score_sets.transfer();
  let ref_2_struct_feature_score_sets = &struct_feature_score_sets;
  let mut bpp_mat_strs = vec![String::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (bpp_mat_str, fasta_record) in multizip((bpp_mat_strs.iter_mut(), fasta_records.iter())) {
      scope.execute(move || {
        let seq_len = fasta_record.seq.len();
        if seq_len <= u8::MAX as usize {
          let bpp_mat = mccaskill_algo::<u8>(
            &fasta_record.seq[..],
            use_contra_model,
            false,
            ref_2_struct_feature_score_sets,
          )
          .0;
          *bpp_mat_str = convert_bpp_mat_2_str(&bpp_mat);
        } else {
          let bpp_mat = mccaskill_algo::<u16>(
            &fasta_record.seq[..],
            use_contra_model,
            false,
            ref_2_struct_feature_score_sets,
          )
          .0;
          *bpp_mat_str = convert_bpp_mat_2_str(&bpp_mat);
        }
      });
    }
  });
  let mut buf_4_writer_2_bpp_mat_file = "# Format = >{RNA sequence id} {line break} {basepairing left nucleotide}, {basepairing right nucleotide}, {basepairing probability} ...".to_string();
  let mut writer_2_bpp_mat_file = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id, bpp_mat_str) in bpp_mat_strs.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    buf_4_rna_id.push_str(bpp_mat_str);
    buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
}

fn convert_bpp_mat_2_str<T>(bpp_mat: &SparseProbMat<T>) -> String
where
  T: HashIndex,
{
  let mut bpp_mat_str = String::new();
  for (&(i, j), &bpp) in bpp_mat.iter() {
    bpp_mat_str.push_str(&format!("{},{},{} ", &i.to_string(), &j.to_string(), bpp));
  }
  bpp_mat_str
}
