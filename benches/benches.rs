extern crate criterion;
extern crate rna_algos;

use criterion::{criterion_group, criterion_main, Criterion};
use rna_algos::durbin_algo::*;
use rna_algos::mccaskill_algo::*;
use rna_algos::utils::*;

fn bench_mccaskill(criterion: &mut Criterion) {
  let fasta_file_reader = Reader::from_file(Path::new(&EXAMPLE_FASTA_FILE_PATH)).unwrap();
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
  let mut struct_feature_score_sets = StructFeatureCountSets::new(0.);
  struct_feature_score_sets.transfer();
  let num_of_threads = num_cpus::get() as NumOfThreads;
  let mut thread_pool = Pool::new(num_of_threads);
  let ref_2_struct_feature_score_sets = &struct_feature_score_sets;
  criterion.bench_function("mccaskill_algo::<u8> (use_contra_model = false)", |b| {
    b.iter(|| {
      thread_pool.scoped(|scope| {
        for fasta_record in &fasta_records {
          scope.execute(move || {
            let _ = mccaskill_algo::<u8>(
              &fasta_record.seq[..],
              false,
              false,
              ref_2_struct_feature_score_sets,
            );
          });
        }
      });
    });
  });
  criterion.bench_function("mccaskill_algo::<u8> (use_contra_model = true)", |b| {
    b.iter(|| {
      thread_pool.scoped(|scope| {
        for fasta_record in &fasta_records {
          scope.execute(move || {
            let _ = mccaskill_algo::<u8>(
              &fasta_record.seq[..],
              true,
              false,
              ref_2_struct_feature_score_sets,
            );
          });
        }
      });
    });
  });
}

fn bench_durbin(criterion: &mut Criterion) {
  let fasta_file_reader = Reader::from_file(Path::new(&EXAMPLE_FASTA_FILE_PATH)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = convert(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    let seq_len = seq.len();
    if seq_len > max_seq_len {
      max_seq_len = seq_len;
    }
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut align_feature_score_sets = AlignFeatureCountSets::new(0.);
  align_feature_score_sets.transfer();
  // let seqs = fasta_records.iter().map(|x| &x.seq[..]).collect();
  let num_of_threads = num_cpus::get() as NumOfThreads;
  let mut thread_pool = Pool::new(num_of_threads);
  let ref_2_align_feature_score_sets = &align_feature_score_sets;
  let num_of_recs = fasta_records.len();
  criterion.bench_function("durbin_algo::<u8>", |b| {
    b.iter(|| {
      thread_pool.scoped(|scope| {
        for i in 0..num_of_recs {
          for j in i + 1..num_of_recs {
            let seq_pair = (&fasta_records[i].seq[..], &fasta_records[j].seq[..]);
            scope.execute(move || {
              let _ = durbin_algo(&seq_pair, ref_2_align_feature_score_sets);
            });
          }
        }
      });
    });
  });
}

criterion_group!(benches, bench_mccaskill, bench_durbin);
criterion_main!(benches);
