extern crate rna_algos;

use rna_algos::mccaskill_algo::*;
use rna_algos::centroid_estimator::*;
use rna_algos::utils::*;

pub type SparseProbMats<T> = Vec<SparseProbMat<T>>;

const MIN_POW_OF_2: i32 = -7;
const MAX_POW_OF_2: i32 = 10;
const DEFAULT_GAMMA: Prob = NEG_INFINITY;

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "An input FASTA file path containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "An output directory path", "STR");
  opts.optopt("g", "gamma", "A specific gamma parameter rather than a range of gamma parameters", "FLOAT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
  opts.optflag("c", "uses_contra_model", "Use the CONTRAfold model instead of Turner's model to score RNA secondary structures");
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
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
  let gamma = if matches.opt_present("gamma") {
    matches.opt_str("gamma").unwrap().parse().unwrap()
  } else {
    DEFAULT_GAMMA
  };
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let uses_contra_model = matches.opt_present("c");
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let seq = convert(fasta_record.seq());
    let seq_len = seq.len();
    if seq_len > max_seq_len {
      max_seq_len = seq_len;
    }
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut thread_pool = Pool::new(num_of_threads);
  if max_seq_len <= u8::MAX as usize {
    multi_threaded_centroid_estimator::<u8>(&mut thread_pool, &fasta_records, output_dir_path, uses_contra_model, gamma);
  } else {
    multi_threaded_centroid_estimator::<u16>(&mut thread_pool, &fasta_records, output_dir_path, uses_contra_model, gamma);
  }
}

fn multi_threaded_centroid_estimator<T>(thread_pool: &mut Pool, fasta_records: &FastaRecords, output_dir_path: &Path, uses_contra_model: bool, gamma: Prob)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let num_of_fasta_records = fasta_records.len();
  let mut bpp_mats = vec![SparseProbMat::<T>::default(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (fasta_record, bpp_mat) in fasta_records.iter().zip(bpp_mats.iter_mut()) {
      let ref seq = fasta_record.seq;
      scope.execute(move || {
        *bpp_mat = mccaskill_algo::<T>(seq, uses_contra_model).0;
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  if gamma != NEG_INFINITY {
    let output_file_path = output_dir_path.join(&format!("gamma={}.fa", gamma));
    compute_and_write_mea_sss::<T>(&bpp_mats, fasta_records, gamma, &output_file_path);
  } else {
    thread_pool.scoped(|scope| {
      for pow_of_2 in MIN_POW_OF_2 .. MAX_POW_OF_2 + 1 {
        let gamma = (2. as Prob).powi(pow_of_2);
        let ref ref_2_bpp_mats = bpp_mats;
        let ref ref_2_fasta_records = fasta_records;
        let output_file_path = output_dir_path.join(&format!("gamma={}.fa", gamma));
        scope.execute(move || {
          compute_and_write_mea_sss::<T>(ref_2_bpp_mats, ref_2_fasta_records, gamma, &output_file_path);
        });
      }
    });
  }
}

fn compute_and_write_mea_sss<T>(bpp_mats: &SparseProbMats<T>, fasta_records: &FastaRecords, gamma: Prob, output_file_path: &Path)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display,
{
  let num_of_fasta_records = fasta_records.len();
  let mut buf = String::new();
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id, fasta_record) in fasta_records.iter().enumerate() {
    let ref bpp_mat = bpp_mats[rna_id];
    let mea_ss = centroid_estimator::<T>(&bpp_mat, fasta_record.seq.len(), gamma);
    let buf_4_rna_id = format!(">{}\n", rna_id) + &unsafe {String::from_utf8_unchecked(get_mea_ss_str::<T>(&mea_ss, fasta_records[rna_id].seq.len()))} + if rna_id < num_of_fasta_records - 1 {"\n"} else {""};
    buf.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_output_file.write_all(buf.as_bytes());
}

fn get_mea_ss_str<T>(mea_ss: &MeaSs<T>, seq_len: usize) -> MeaSsStr
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let mut mea_ss_str = vec![UNPAIRING_BASE; seq_len];
  for &(i, j) in &mea_ss.bp_pos_pairs {
    mea_ss_str[i.to_usize().unwrap()] = BASE_PAIRING_LEFT_BASE;
    mea_ss_str[j.to_usize().unwrap()] = BASE_PAIRING_RIGHT_BASE;
  }
  mea_ss_str
}
