extern crate criterion;
extern crate rna_algos;

use criterion::{criterion_group, criterion_main, Criterion};
use rna_algos::durbin_algo::*;
use rna_algos::mccaskill_algo::*;
use rna_algos::utils::*;

fn bench_mccaskill_algo(criterion: &mut Criterion) {
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
  let uses_contra_model = false;
  let allows_short_hairpins = false;
  criterion.bench_function("mccaskill_algo_algo::<u8> (uses_contra_model = false)", |x| {
    let y = &fold_scores;
    x.iter(|| {
      thread_pool.scoped(|z| {
        for a in &fasta_records {
          z.execute(move || {
            let _ = mccaskill_algo::<u8>(
              &a.seq[..],
              uses_contra_model,
              allows_short_hairpins,
              y,
            );
          });
        }
      });
    });
  });
  let uses_contra_model = true;
  criterion.bench_function("mccaskill_algo::<u8> (uses_contra_model = true)", |x| {
    let y = &fold_scores;
    x.iter(|| {
      thread_pool.scoped(|z| {
        for a in &fasta_records {
          z.execute(move || {
            let _ = mccaskill_algo::<u8>(
              &a.seq[..],
              uses_contra_model,
              allows_short_hairpins,
              y,
            );
          });
        }
      });
    });
  });
}

fn bench_durbin_algo(criterion: &mut Criterion) {
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
  criterion.bench_function("durbin_algo::<u8>", |x| {
    let y = &align_scores;
    x.iter(|| {
      thread_pool.scoped(|z| {
        for i in 0..num_records {
          for j in i + 1..num_records {
            let a = (&fasta_records[i].seq[..], &fasta_records[j].seq[..]);
            z.execute(move || {
              let _ = durbin_algo(&a, y);
            });
          }
        }
      });
    });
  });
}

criterion_group!(benches, bench_mccaskill_algo, bench_durbin_algo);
criterion_main!(benches);
