#[allow(unused_imports)]
use rustc_hash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use crossbeam::queue::ArrayQueue;
use seq_io::fasta::{Reader, Record};
use bio::alignment::pairwise::*;
use wyhash::WyHash;
use std::hash::{Hash, Hasher};

pub mod utils;
use utils::align::aligner;

static BYTE2BITS: [u8; 256] = {
    let mut l = [0u8; 256];
    l[b'A' as usize] = 0b00;
    l[b'a' as usize] = 0b00;
    l[b'C' as usize] = 0b01;
    l[b'c' as usize] = 0b01;
    l[b'G' as usize] = 0b10;
    l[b'g' as usize] = 0b10;
    l[b'T' as usize] = 0b11;
    l[b't' as usize] = 0b11;
    l
};

fn hash(seed: u64, kmer: u128) -> u64 {
    let mut h = WyHash::with_seed(seed);
    kmer.hash(&mut h);
    h.finish()
}

fn sketch<T: AsRef<[u8]>>(seq: &T, k: usize) -> FxHashMap<u64,Vec<usize>> {
    let mut sketch: FxHashMap<u64,Vec<usize>> = FxHashMap::default();
    let mask : u128 = (1 << 2*k) - 1;
    let seq = seq.as_ref();
    for i in 0..seq.len()-k {
        let kmer = &seq[i..i+k];
        let mut kmer_b: u128 = 0;
        for b in kmer {
            unsafe {
                kmer_b = (kmer_b << 2) & mask | ((*BYTE2BITS.as_ptr().add(*b as usize)) as u128);
            }
        }
        let hash = hash(0, kmer_b);
        sketch.entry(hash).or_insert(Vec::new()).push(i);
    }
    return sketch;
}

#[allow(dead_code)]
fn main() {
    let _k: usize = 15;
    let _w: usize = 138;
    let ref_path = "/usr1/aidanz/projects/read_sim/data/fasta/U00096.3.fasta";
    let mut ref_reader = Reader::from_path(ref_path).unwrap();
    let ref_record = ref_reader.next().expect("Failed to read reference record").unwrap();
    let ref_seq = ref_record.owned_seq();
    let ref_sketch = sketch(&ref_seq, _k);


    let probe_path = "/usr1/aidanz/projects/read_sim/data/fasta/all_probes.fa";
    let mut probe_reader = Reader::from_path(probe_path).unwrap();
    let mut probe_sketches : Vec<(Vec<u8>,FxHashMap<u64,Vec<usize>>)> = Vec::new();

    println!("Reading probes");
    while let Some(record) = probe_reader.next() {
        let record = record.expect("Error reading record");
        let seq = record.owned_seq();
        probe_sketches.push((seq.clone(),sketch(&seq, _k)));
    }
    println!("Finished reading probes");

    let score = |a: u8, b: u8| if a == b { 2i32 } else { -2i32 };

    let mut pileup : FxHashMap<usize,i32> = FxHashMap::default();

    let q : ArrayQueue<(usize,usize,usize)> = ArrayQueue::new(probe_sketches.len());

    (0..probe_sketches.len()).into_par_iter().for_each(|i| {
        println!("Processing probe {}", i);
        let mut aligner = Aligner::new(-4,-1,&score);
        //let mut aligner = aligner::new(0,0,-3,-1,2,-2);
        let sk = &probe_sketches[i].1;
        let seq = &probe_sketches[i].0;
        for hash_key in sk.keys() {
            match ref_sketch.get(hash_key) {
                Some(ref_positions) => {
                    for ref_pos in ref_positions {
                        let ref_start = ref_pos.saturating_sub(_w);
                        let ref_end = (*ref_pos + _w).min(ref_seq.len());
                        let ref_window = &ref_seq[ref_start..ref_end];                         
                        let alignment = aligner.semiglobal(&seq, &ref_window);
                        if alignment.score >= 80 {
                            let _ = q.push((i,
                                    alignment.ystart as usize + ref_start,
                                    alignment.ylen as usize));
                        }
                    }
                },
                None => {}
            }
        }
    });

    let probe_hits : FxHashSet<(usize,usize,usize)> = FxHashSet::from_iter(q.into_iter());

    for hit in probe_hits.iter() {
        println!("Probe {} hit at position {} with length {}", hit.0, hit.1, hit.2);
        for i in hit.1..hit.2+hit.1 {
            *pileup.entry(i).or_insert(0) += 1;
        }
    }
    let pileup_count = pileup.len();
    println!("Number of positions with at least one hit: {}", pileup_count);
    println!("Number of probes with at least one hit: {}", probe_hits.len());
}
