use rustc_hash::{FxHashMap, FxHashSet};
use seq_io::fasta::{Reader, Record};
use wyhash::WyHash;
use std::any;
use std::hash::{Hash, Hasher};
use std::fs::File;
use std::io::{BufWriter, BufReader, prelude::*};

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

fn sketch<T: AsRef<[u8]>>(seq: T, k: usize) -> FxHashMap<u64,Vec<usize>> {
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

fn main() {
    let k: usize = 15;
    let ref_path = "/usr1/aidanz/projects/read_sim/data/fasta/U00096.3.fasta";
    let mut ref_reader = Reader::from_path(ref_path).unwrap();
    let ref_record = ref_reader.next().expect("Failed to read reference record").unwrap();
    let ref_sketch = sketch(ref_record.owned_seq(), k);


    let probe_path = "/usr1/aidanz/projects/read_sim/data/fasta/all_probes.fa";
    let mut probe_reader = Reader::from_path(probe_path).unwrap();
    let mut probe_sketches : Vec<FxHashMap<u64,Vec<usize>>> = Vec::new();

    println!("Reading probes");
    while let Some(record) = probe_reader.next() {
        let record = record.expect("Error reading record");
        let seq = record.owned_seq();
        probe_sketches.push(sketch(seq, k));
    }

    let mut count = 0;
    for sk in probe_sketches.iter() {
        let mut matches = 0;
        for hash_key in sk.keys() {
            if ref_sketch.contains_key(hash_key) {
                matches += 1;
            }
        }
        if matches > 0 {
            println!("Found {} matching k-mers", matches);
            count += 1;
        }
    }
    println!("Found {} probes with matching k-mers", count);
}
