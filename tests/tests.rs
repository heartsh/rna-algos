extern crate rna_algos;

use rna_algos::durbin_algo::*;
use rna_algos::mccaskill_algo::*;
use rna_algos::utils::*;

#[test]
fn test_mccaskill() {
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
  thread_pool.scoped(|scope| {
    for fasta_record in &fasta_records {
      scope.execute(move || {
        let bpp_mat = mccaskill_algo::<u8>(
          &fasta_record.seq[..],
          false,
          false,
          ref_2_struct_feature_score_sets,
        )
        .0;
        for &bpp in bpp_mat.values() {
          assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&bpp));
        }
        let bpp_mat = mccaskill_algo::<u8>(
          &fasta_record.seq[..],
          true,
          false,
          ref_2_struct_feature_score_sets,
        )
        .0;
        for &bpp in bpp_mat.values() {
          assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&bpp));
        }
      });
    }
  });
}

#[test]
fn test_durbin() {
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
  let num_of_threads = num_cpus::get() as NumOfThreads;
  let mut thread_pool = Pool::new(num_of_threads);
  let ref_2_align_feature_score_sets = &align_feature_score_sets;
  let num_of_recs = fasta_records.len();
  thread_pool.scoped(|scope| {
    for i in 0..num_of_recs {
      for j in i + 1..num_of_recs {
        let seq_pair = (&fasta_records[i].seq[..], &fasta_records[j].seq[..]);
        scope.execute(move || {
          let align_prob_mat = durbin_algo(&seq_pair, ref_2_align_feature_score_sets);
          for &align_prob in align_prob_mat.iter().flatten() {
            assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&align_prob));
          }
        });
      }
    }
  });
}
