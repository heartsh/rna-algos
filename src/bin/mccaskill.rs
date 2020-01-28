extern crate rna_algos;
extern crate scoped_threadpool;
extern crate bio;
extern crate num_cpus;
extern crate itertools;

use rna_algos::mccaskill_algo::*;
use rna_algos::utils::*;
use scoped_threadpool::Pool;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::fs::create_dir;

type NumOfThreads = u32;

const BPP_MAT_FILE_NAME: &'static str = "bpp_mats.dat";
const UPP_MAT_FILE_NAME: &'static str = "upp_mats.dat";
const VERSION: &'static str = "0.1.8";

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_file_path = opts.opt_str("i").expect("Failed to get the path to an input FASTA file containing RNA sequences from command arguments.");
  let input_file_path = Path::new(&input_file_path);
  let output_dir_path = opts.opt_str("o").expect("Failed to get the path to an output directory from command arguments.");
  let output_dir_path = Path::new(&output_dir_path);
  let num_of_threads = if opts.opt_present("t") {
    opts.opt_str("t").expect("Failed to get the number of threads in multithreading from command arguments.").parse().expect("Failed to parse the number of threads in multithreading.")
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).expect("Failed to set a FASTA file reader.");
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.expect("Failed to read a FASTA record.");
    let seq = unsafe {from_utf8_unchecked(fasta_record.seq()).to_uppercase().as_bytes().iter().filter(|&&base| {is_rna_base(base)}).map(|&base| {base}).collect()};
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let num_of_fasta_records = fasta_records.len();
  let mut thread_pool = Pool::new(num_of_threads);
  let mut bpp_mats = vec![SparseProbMat::default(); num_of_fasta_records];
  let mut upp_mats = vec![Probs::new(); num_of_fasta_records];
  let mut max_free_energies = vec![0.; num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (bpp_mat, upp_mat, max_free_energy, fasta_record) in multizip((bpp_mats.iter_mut(), upp_mats.iter_mut(), max_free_energies.iter_mut(), fasta_records.iter())) {
      scope.execute(move || {
        let (obtained_bpp_mat, obtained_upp_mat, obtained_max_free_energy) = get_bpp_and_unpair_prob_mats(&fasta_record.seq[..]);
        *bpp_mat = obtained_bpp_mat;
        *upp_mat = obtained_upp_mat;
        *max_free_energy = obtained_max_free_energy;
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let mut buf_4_writer_2_bpp_mat_file = format!("; The version {} of the McCaskill program.\n; The path to the input file in order to compute the base-pairing Probability matrices on secondary structure in this file = \"{}\".\n; The values of the parameters used in order to compute the matrices are as follows.\n; \"num_of_threads\" = {}.", VERSION, input_file_path.display(), num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence and (the scaled inverse of) the minimum free energy of the sequence on secondary structure. The next row to the row is with the base-pairing probability matrix on secondary structure of the sequence.";
  let bpp_mat_file_path = output_dir_path.join(BPP_MAT_FILE_NAME);
  let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).expect("Failed to create an output file."));
  for (rna_id, bpp_mat) in bpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in bpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i, j, bpp));
    }
    buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
  let mut buf_4_writer_2_upp_mat_file = format!("; The version {} of the McCaskill program.\n; The path to the input file in order to compute the UnPairing probability matrices on secondary structure in this file = \"{}\".\n; The values of the parameters used in order to compute the matrices are as follows.\n; \"num_of_threads\" = {}.", VERSION, input_file_path.display(), num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The next row to the row is with the unpairing probability matrix on secondary structure of the sequence on secondary structure.";
  let upp_mat_file_path = output_dir_path.join(UPP_MAT_FILE_NAME);
  let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).expect("Failed to create an output file."));
  for (rna_id, upp_mat) in upp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (i, &upp) in upp_mat.iter().enumerate() {
      buf_4_rna_id.push_str(&format!("{},{} ", i, upp));
    }
    buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
}
