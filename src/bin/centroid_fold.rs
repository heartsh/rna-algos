extern crate rna_algos;

use rna_algos::centroid_fold::*;
use rna_algos::mccaskill_algo::*;
use rna_algos::utils::*;

pub type SparseProbMats<T> = Vec<SparseProbMat<T>>;

const MIN_POW_2: i32 = -7;
const MAX_POW_2: i32 = 10;
const DEFAULT_CENTROID_THRESHOLD: Prob = NEG_INFINITY;

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
  opts.reqopt("o", "output_dir_path", "An output directory path", "STR");
  opts.optopt(
    "g",
    "centroid_threshold",
    "A specific centroid threshold rather than a range of centroid thresholds",
    "FLOAT",
  );
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
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
  let centroid_threshold = if matches.opt_present("centroid_threshold") {
    matches.opt_str("centroid_threshold").unwrap().parse().unwrap()
  } else {
    DEFAULT_CENTROID_THRESHOLD
  };
  let num_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumThreads
  };
  let uses_contra_model = matches.opt_present("c");
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let seq = bytes2seq(fasta_record.seq());
    let seq_len = seq.len();
    if seq_len > max_seq_len {
      max_seq_len = seq_len;
    }
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut thread_pool = Pool::new(num_threads);
  if max_seq_len <= u8::MAX as usize {
    multithread_centroid_fold::<u8>(
      &mut thread_pool,
      &fasta_records,
      output_dir_path,
      uses_contra_model,
      centroid_threshold,
    );
  } else {
    multithread_centroid_fold::<u16>(
      &mut thread_pool,
      &fasta_records,
      output_dir_path,
      uses_contra_model,
      centroid_threshold,
    );
  }
}

fn multithread_centroid_fold<T>(
  thread_pool: &mut Pool,
  fasta_records: &FastaRecords,
  output_dir_path: &Path,
  uses_contra_model: bool,
  centroid_threshold: Prob,
) where
  T: HashIndex,
{
  let mut fold_score_sets = FoldScoreSets::new(0.);
  fold_score_sets.transfer();
  let ref_fold_score_sets = &fold_score_sets;
  let num_fasta_records = fasta_records.len();
  let mut basepair_prob_mats = vec![SparseProbMat::<T>::default(); num_fasta_records];
  let allows_short_hairpins = false;
  thread_pool.scoped(|scope| {
    for (fasta_record, basepair_probs) in fasta_records.iter().zip(basepair_prob_mats.iter_mut()) {
      let seq = &fasta_record.seq;
      scope.execute(move || {
        *basepair_probs = mccaskill_algo::<T>(
          seq,
          uses_contra_model,
          allows_short_hairpins,
          ref_fold_score_sets,
        )
        .0;
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  if centroid_threshold != NEG_INFINITY {
    let output_file_path = output_dir_path.join(format!("centroid_threshold={}.fa", centroid_threshold));
    write_centroid_fold::<T>(&basepair_prob_mats, fasta_records, centroid_threshold, &output_file_path);
  } else {
    let ref_basepair_prob_mats = &basepair_prob_mats;
    thread_pool.scoped(|scope| {
      for pow_2 in MIN_POW_2..MAX_POW_2 + 1 {
        let centroid_threshold  = (2. as Prob).powi(pow_2);
        let output_file_path = output_dir_path.join(format!("centroid_threshold={}.fa", centroid_threshold));
        scope.execute(move || {
          write_centroid_fold::<T>(ref_basepair_prob_mats, fasta_records, centroid_threshold, &output_file_path);
        });
      }
    });
  }
}

fn write_centroid_fold<T>(
  basepair_prob_mats: &SparseProbMats<T>,
  fasta_records: &FastaRecords,
  centroid_threshold: Prob,
  output_file_path: &Path,
) where
  T: HashIndex,
{
  let num_fasta_records = fasta_records.len();
  let mut buf = String::new();
  let mut writer = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id, fasta_record) in fasta_records.iter().enumerate() {
    let basepair_probs = &basepair_prob_mats[rna_id];
    let centroid_fold = centroid_fold::<T>(basepair_probs, fasta_record.seq.len(), centroid_threshold);
    let buf_rna_id = format!(">{}\n", rna_id)
      + &unsafe {
        String::from_utf8_unchecked(get_fold_str::<T>(
          &centroid_fold.basepair_pos_pairs,
          fasta_records[rna_id].seq.len(),
        ))
      }
      + if rna_id < num_fasta_records - 1 {
        "\n"
      } else {
        ""
      };
    buf.push_str(&buf_rna_id);
  }
  let _ = writer.write_all(buf.as_bytes());
}

fn get_fold_str<T>(fold: &PosPairs<T>, seq_len: usize) -> FoldStr
where
  T: HashIndex,
{
  let mut fold_str = vec![UNPAIR; seq_len];
  for &(i, j) in fold {
    fold_str[i.to_usize().unwrap()] = BASEPAIR_LEFT;
    fold_str[j.to_usize().unwrap()] = BASEPAIR_RIGHT;
  }
  fold_str
}
