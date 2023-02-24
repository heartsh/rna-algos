extern crate rna_algos;

use rna_algos::durbin_algo::*;
use rna_algos::mccaskill_algo::*;
use rna_algos::utils::*;

#[test]
fn test_mccaskill_algo() {
  let fasta_file_reader = Reader::from_file(Path::new(&EXAMPLE_FASTA_FILE_PATH)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for x in fasta_file_reader.records() {
    let x = x.unwrap();
    let y = bytes2seq(x.seq());
    let z = y.len();
    if z > max_seq_len {
      max_seq_len = z;
    }
    fasta_records.push(FastaRecord::new(String::from(x.id()), y));
  }
  let mut fold_scores = FoldScoreSets::new(0.);
  fold_scores.transfer();
  let num_threads = num_cpus::get() as NumThreads;
  let mut thread_pool = Pool::new(num_threads);
  let allows_short_hairpins = false;
  thread_pool.scoped(|x| {
    let y = &fold_scores;
    for z in &fasta_records {
      x.execute(move || {
        let uses_contra_model = false;
        let a = mccaskill_algo::<u8>(&z.seq[..], uses_contra_model, allows_short_hairpins, y).0;
        for &a in a.values() {
          assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&a));
        }
        let uses_contra_model = true;
        let a = mccaskill_algo::<u8>(&z.seq[..], uses_contra_model, allows_short_hairpins, y).0;
        for &a in a.values() {
          assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&a));
        }
      });
    }
  });
}

#[test]
fn test_durbin_algo() {
  let fasta_file_reader = Reader::from_file(Path::new(&EXAMPLE_FASTA_FILE_PATH)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for x in fasta_file_reader.records() {
    let x = x.unwrap();
    let mut y = bytes2seq(x.seq());
    y.insert(0, PSEUDO_BASE);
    y.push(PSEUDO_BASE);
    let z = y.len();
    if z > max_seq_len {
      max_seq_len = z;
    }
    fasta_records.push(FastaRecord::new(String::from(x.id()), y));
  }
  let mut align_scores = AlignScores::new(0.);
  align_scores.transfer();
  let num_threads = num_cpus::get() as NumThreads;
  let mut thread_pool = Pool::new(num_threads);
  let num_records = fasta_records.len();
  thread_pool.scoped(|x| {
    let y = &align_scores;
    for i in 0..num_records {
      for j in i + 1..num_records {
        let z = (&fasta_records[i].seq[..], &fasta_records[j].seq[..]);
        x.execute(move || {
          let x = durbin_algo(&z, y);
          for &x in x.iter().flatten() {
            assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&x));
          }
        });
      }
    }
  });
}
