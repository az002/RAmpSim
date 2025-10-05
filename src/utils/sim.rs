use std::collections::{VecDeque, HashMap};
use bam::{ header, Record as BamRecord};
use rand::prelude::*;
use rand_distr::{LogNormal, Distribution};
use std::fs::File;
use std::io::{BufWriter, BufReader};
use bam::SamReader;
use nalgebra::ArrayStorage;
use nalgebra::SVector;
use seq_io::fasta::{Reader, Record};
use std::cmp;
use std::io::Write;

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

const BETA: f64 = 1.0/0.001987;

#[inline]
fn nt2bits(b: u8) -> u8 { BYTE2BITS[b as usize] }

#[inline]
const fn comp2(x: u8) -> u8 { x ^ 0b11 } // A<->T, C<->G on 2-bit alphabet

#[inline]
const fn revcomp3(k: u8) -> u8 {
    let a = comp2((k >> 4) & 0b11);
    let b = comp2((k >> 2) & 0b11);
    let c = comp2(k & 0b11);
    (c << 4) | (b << 2) | a
}

const fn build_trimer_index() -> [u8; 64] {
	let mut map = [u8::MAX; 64];
    let mut next = 0u8;
    let mut k = 0u8;
    while k < 64 {
        let c = {
            let rc = revcomp3(k);
            if rc < k { rc } else { k }
        };
        if map[c as usize] == u8::MAX {
            map[c as usize] = next;
            next += 1;
        }
        map[k as usize] = map[c as usize];
        k += 1;
    }
    map
}

const TRIMER_INDEX: [u8; 64] = build_trimer_index();

const SCORE : SVector<f64, 32> = SVector::<f64,32>::from_array_storage(ArrayStorage([[-1.0333523772924682,  -1.458111012202453,  -1.035949547554549,  -0.7124736682600892, 
    -1.8355462789826653,  -3.2816907774428596,  -2.9387615700865912,  -2.4322581762379056, 
    -6.050351136676882e-26,  -1.16892836679913,  -1.2733818523786504,  -0.42044691866600825, 
    -2.2817669581121254,  -1.4822900733089985,  -0.9125362005565887,  -3.670822543027401e-23, 
    -1.273598247757921,  -1.3530081563908097e-28,  -0.19137969659843326,  -1.2926131142356785, 
    -0.3257778299362714,  -2.1797369472925356,  -1.3319137918104318,  -1.7767543614895418, 
    -2.617678120625974,  -1.5216318279285466e-24,  -0.8229383115364808,  -2.1347399492841697, 
    -1.502770970281663e-38,  -0.8507439536421594,  -0.31680823162220356,  -0.06611591220578053]]));

pub struct Hit {
	pub score: f64,
	pub seq_id: String,
	pub start: usize,
	pub end: usize,
}

pub struct SequenceInput {
	pub samreader: SamReader<BufReader<File>>,
	pub output_writer: BufWriter<File>,
}


pub struct Simulator {
    fragment_len: usize,
	lognorm_sd: f64,
	nfragments: usize,
	split: f64,
	temperature: f64,
	probe_multiplicities: HashMap<String, usize>,
	abundances: HashMap<String, f64>,
	seqid2refid: HashMap<String, String>,
}

impl Simulator {
	pub fn new(
		fragment_len: usize, 
		lognorm_sd: f64, 
		nfragments: usize,
		split: f64,
		temperature: f64,
		probe_multiplicities: HashMap<String, usize>,
		abundances: HashMap<String, f64>,
		seqid2refid: HashMap<String, String>,
	) -> Self {
		let mut sim = Simulator {
			fragment_len,
			lognorm_sd,
			nfragments,
			split,
			temperature,
			probe_multiplicities,
			abundances,
			seqid2refid,
		};
		sim
	}

	pub fn process_hit(record: &BamRecord, seqid: String, temp: f64, abundances: &HashMap<String, f64>, seqid2refid: &HashMap<String, String>) -> Hit{
		let mut ftrimer = 0u8;
		let mut rtrimer = 0u8;
		let shift = 4;
        let mask:u8 = 0b00111111;
        let mut profile = SVector::<f64, 32>::zeros();
        let mut count : u8 = 0;
        for entry in record.alignment_entries().unwrap() {
            if entry.is_seq_match() {
				ftrimer = (ftrimer << 2) & mask | nt2bits(entry.record_nt().unwrap());
				rtrimer = (rtrimer >> 2) | (comp2(nt2bits(entry.record_nt().unwrap())) << shift);
				let trimer = cmp::min(ftrimer,rtrimer);
                count += 1;
                if count >= 3 {
                    profile[TRIMER_INDEX[trimer as usize] as usize] += 1.;
                }
            }
            else {
                ftrimer = 0;
				rtrimer = 0;
                count = 0;
            }
        }
		// println!("Processing record: {}", seqid);
		let refname = seqid2refid.get(&seqid).expect("Reference name not found");
		let abundance = *abundances.get(refname).unwrap_or(&0.0);
        Hit {
            score: (-BETA/temp*profile.dot(&SCORE)).exp() * abundance,
            seq_id: seqid,
            start: record.start() as usize,
            end: record.calculate_end() as usize,
        }
	}

	pub fn compute_telseq_distr(&mut self, seqinput: &mut SequenceInput) -> HashMap<String, Vec<Hit>> {
		let mut probe_hits: HashMap<String, Vec<Hit>> = HashMap::new();
		let header = seqinput.samreader.header().to_owned();
		for result in &mut seqinput.samreader {
			match result {
				Ok(record) => {
					let seqid = header.reference_name(record.ref_id() as u32).unwrap().to_string();
					let hit = Simulator::process_hit(&record, seqid, self.temperature, &self.abundances, &self.seqid2refid);
					probe_hits.entry(String::from_utf8(record.name().to_vec()).expect("Invalid UTF-8")).or_default().push(hit);
				},
				Err(e) => {
					eprintln!("Error reading SAM record: {}", e);
					continue;
				}
			}
		}
		probe_hits
	}

	pub fn sample_telseq(&mut self, seqinput: &mut SequenceInput, ref_seqs : &HashMap<String, String>) {
		let mut probe_hits = self.compute_telseq_distr(seqinput);
		let probs: Vec<f64> = self.probe_multiplicities.iter()
			.filter(|(k, _)| probe_hits.contains_key(*k))
			.map(|(_, v)| *v as f64)
			.collect();
		let probe_names: Vec<String> = self.probe_multiplicities.keys()
			.filter(|k| probe_hits.contains_key(*k))
			.cloned()
			.collect();

		let mut rng = rand::rng();
		let multinomial = rand::distr::weighted::WeightedIndex::new(&probs).unwrap();

		let lognorm = LogNormal::new((self.fragment_len as f64).ln(), self.lognorm_sd).unwrap();
		let n = (self.nfragments as f64 * self.split) as usize;
		eprintln!("Sampling {} telseq fragments", n);
		for _ in 0..n {
			let idx = multinomial.sample(&mut rng);
			let probe_name = &probe_names[idx];

			if let Some(hits) = probe_hits.get_mut(probe_name) {
				if !hits.is_empty() {
					let hit_probs = hits.iter().map(|h| h.score).collect::<Vec<f64>>();
					let hit_dist = rand::distr::weighted::WeightedIndex::new(&hit_probs).unwrap();
					let hit_idx = hit_dist.sample(&mut rng);
					let hit = &hits[hit_idx];
					let ref_seq = ref_seqs.get(&hit.seq_id).expect("Reference sequence not found");

					let frag_len = lognorm.sample(&mut rng) as usize;

					let min_center = hit.start.saturating_sub(self.fragment_len);
					let max_center = hit.end + self.fragment_len;
					let max_center = cmp::min(max_center, ref_seq.len());

					let center = if max_center > min_center {
						rng.random_range(min_center..max_center)
					} else {
						min_center
					};

					let frag_start = center.saturating_sub(frag_len / 2);
					let frag_end = cmp::min(frag_start + frag_len, ref_seq.len());

					// Write to output
					let fragment = &ref_seq[frag_start..frag_end];
					writeln!(seqinput.output_writer, ">{}|{}|{}|{}-{}\n{}",
						probe_name, hit.seq_id, hit.start, frag_start, frag_end, fragment
					).expect("Error writing to output");

					println!("[TELSEQ]:{}:{}:{}-{}:{}",self.seqid2refid.get(&hit.seq_id).unwrap(), hit.seq_id, frag_start, frag_end, probe_name);
				}
			}
		}
	}

}