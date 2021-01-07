extern crate rna_algos;

pub use rna_algos::mccaskill_algo::*;
pub use rna_algos::utils::*;


fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_file_path", "The path to an output file", "STR");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
  opts.optflag("c", "uses_contra_model", "Use CONTRAfold model instead of Turner's model to score RNA secondary structures");
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
  let uses_contra_model = matches.opt_present("c");
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let seq = convert(fasta_record.seq());
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let num_of_fasta_records = fasta_records.len();
  let mut thread_pool = Pool::new(num_of_threads);
  let mut bpp_mat_strs = vec![String::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (bpp_mat_str, fasta_record) in multizip((bpp_mat_strs.iter_mut(), fasta_records.iter())) {
      scope.execute(move || {
        let seq_len = fasta_record.seq.len();
        if seq_len <= u8::MAX as usize {
          let bpp_mat = mccaskill_algo::<u8>(&fasta_record.seq[..], uses_contra_model);
          *bpp_mat_str = convert_bpp_mat_2_str(&bpp_mat);
        } else {
          let bpp_mat = mccaskill_algo::<u16>(&fasta_record.seq[..], uses_contra_model);
          *bpp_mat_str = convert_bpp_mat_2_str(&bpp_mat);
        }
      });
    }
  });
  let mut buf_4_writer_2_bpp_mat_file = format!("; The path to the input file in order to compute the base-pairing Probability matrices on secondary structure in this file = \"{}\".\n; The values of the parameters used in order to compute the matrices are as follows.\n; \"num_of_threads\" = {}.", input_file_path.display(), num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The next row to the row is with the base-pairing probability matrix on secondary structure of the sequence.";
  let mut writer_2_bpp_mat_file = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id, bpp_mat_str) in bpp_mat_strs.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    buf_4_rna_id.push_str(&bpp_mat_str);
    buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
}

fn convert_bpp_mat_2_str<T>(bpp_mat: &SparseProbMat<T>) -> String
where
    T: Display + Copy,
{
  let mut bpp_mat_str = String::new();
  for (&(i, j), &bpp) in bpp_mat.iter() {
    bpp_mat_str.push_str(&format!("{},{},{} ", &i.to_string(), &j.to_string(), bpp));
  }
  bpp_mat_str
}
