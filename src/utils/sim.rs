use std::collections::HashMap;
use bam::Record as BamRecord;
use rand::prelude::*;
use rand_distr::{LogNormal, Distribution};
use std::fs::File;
use std::io::{BufWriter, BufReader};
use bam::SamReader;
use std::cmp;
use std::io::Write;
use crate::utils::thermo;

pub struct Hit {
	pub score: f64,
	pub seq_id: String,
	pub start: usize,
	pub end: usize,
}


pub struct Simulator {
    fragment_len: usize,
	lognorm_sd: f64,
	nfragments: usize,
	sampled: usize,
	split: f64,
	temperature: f64,
	probe_multiplicities: HashMap<String, usize>,
	abundances: HashMap<String, f64>,
	seqid2refid: HashMap<String, String>,
	hits: HashMap<String, Vec<Hit>>,
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
		let sim = Simulator {
			fragment_len,
			lognorm_sd,
			nfragments,
			sampled: 0,
			split,
			temperature,
			probe_multiplicities,
			abundances,
			seqid2refid,
			hits: HashMap::new(),
		};
		sim
	}

	pub fn process_hit(record: &BamRecord, seqid: String, temp: f64, abundances: &HashMap<String, f64>, seqid2refid: &HashMap<String, String>) -> Hit{
		let score = thermo::score_bam_alignment_multiT(record, temp, thermo::MismatchStrategy::SkipStacking);
		let refname = seqid2refid.get(&seqid).expect("Reference name not found");
		let abundance = *abundances.get(refname).unwrap_or(&0.0);
        Hit {
            score: score * abundance,
            seq_id: seqid,
            start: record.start() as usize,
            end: record.calculate_end() as usize,
        }
	}

	pub fn compute_telseq_distr(&mut self, seqinput: &mut SamReader<BufReader<File>>) {
		let header = seqinput.header().to_owned();
		for result in seqinput {
			match result {
				Ok(record) => {
					let seqid = header.reference_name(record.ref_id() as u32).unwrap().to_string();
					let hit = Simulator::process_hit(&record, seqid, self.temperature, &self.abundances, &self.seqid2refid);
					self.hits.entry(String::from_utf8(record.name().to_vec()).expect("Invalid UTF-8")).or_default().push(hit);
				},
				Err(e) => {
					eprintln!("Error reading SAM record: {}", e);
					continue;
				}
			}
		}
	}

	pub fn sample_telseq(&mut self, seqinput: &mut BufWriter<File>, ref_seqs : &HashMap<String, String>) {
		// If probe_multiplicities is empty, use uniform distribution over all probes with hits
		let (probs, probe_names): (Vec<f64>, Vec<String>) = if self.probe_multiplicities.is_empty() {
			let names: Vec<String> = self.hits.keys().cloned().collect();
			let uniform_probs: Vec<f64> = vec![1.0; names.len()];
			(uniform_probs, names)
		} else {
			let probs: Vec<f64> = self.probe_multiplicities.iter()
				.filter(|(k, _)| self.hits.contains_key(*k))
				.map(|(_, v)| *v as f64)
				.collect();
			let names: Vec<String> = self.probe_multiplicities.keys()
				.filter(|k| self.hits.contains_key(*k))
				.cloned()
				.collect();
			(probs, names)
		};

		if probs.is_empty() {
			eprintln!("Warning: No probes with hits found. Skipping telseq sampling.");
			return;
		}

		let mut rng = rand::rng();
		let multinomial = rand::distr::weighted::WeightedIndex::new(&probs).unwrap();

		let lognorm = LogNormal::new((self.fragment_len as f64).ln(), self.lognorm_sd).unwrap();
		let n = (self.nfragments as f64 * self.split) as usize;
		eprintln!("Sampling {} telseq fragments", n);
		for _ in 0..n {
			let idx = multinomial.sample(&mut rng);
			let probe_name = &probe_names[idx];

			if let Some(hits) = self.hits.get_mut(probe_name) {
				if !hits.is_empty() {
					let hit_probs = hits.iter().map(|h| h.score).collect::<Vec<f64>>();
					let hit_dist = rand::distr::weighted::WeightedIndex::new(&hit_probs).unwrap();
					let hit_idx = hit_dist.sample(&mut rng);
					let hit = &hits[hit_idx];
					let ref_seq = ref_seqs.get(&hit.seq_id).expect("Reference sequence not found");

					let frag_len = lognorm.sample(&mut rng) as usize;

					let min_center = hit.end.saturating_sub(self.fragment_len/2);
					let max_center = hit.start + self.fragment_len/2;
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
					writeln!(seqinput, ">{}|{}|{}|{}-{}\n{}",
						probe_name, hit.seq_id, hit.start, frag_start, frag_end, fragment
					).expect("Error writing to output");

					println!("[TELSEQ]:{}:{}:{}-{}:{}",self.seqid2refid.get(&hit.seq_id).unwrap(), hit.seq_id, frag_start, frag_end, probe_name);
					self.sampled += 1;
				}
			}
		}
	}
	

	pub fn sample_background(
		&mut self,
		seq_input: &mut BufWriter<File>,
		ref_seqs : &HashMap<String, String>
	) {
		let n_remain = self.nfragments - self.sampled;
		eprintln!("Sampling {} background fragments", n_remain);

		let mut seq_list: Vec<&String> = Vec::new();
		let mut weights: Vec<f64> = Vec::new();

		for (seq_id, ref_id) in &self.seqid2refid {
			if let Some(seq) = ref_seqs.get(seq_id) {
				let abundance = *self.abundances.get(ref_id).unwrap_or(&0.0);
				let weight = abundance * seq.len() as f64;
				if weight > 0.0 {
					seq_list.push(seq_id);
					weights.push(weight);
				}
			}
		}

		if weights.is_empty() {
			eprintln!("Warning: No valid sequences for background sampling.");
			return;
		}

		let dist = rand::distr::weighted::WeightedIndex::new(&weights).unwrap();
		let lognorm = LogNormal::new((self.fragment_len as f64).ln(), self.lognorm_sd).unwrap();
		let mut rng = rand::rng();

		for _ in 0..n_remain {
			let idx = dist.sample(&mut rng);
			let seq_id = seq_list[idx];
			let seq = ref_seqs.get(seq_id).unwrap();

			let frag_len = cmp::min(lognorm.sample(&mut rng) as usize, seq.len());
			let frag_start = rng.random_range(0..=seq.len() - frag_len);
			let frag_end = frag_start + frag_len;

			writeln!(
				seq_input,
				">Background|{}|{}|{}-{}\n{}",
				self.seqid2refid.get(seq_id).unwrap(),
				seq_id,
				frag_start,
				frag_end,
				&seq[frag_start..frag_end]
			).unwrap();

			println!("[BG]:{}:{}:{}-{}:", self.seqid2refid.get(seq_id).unwrap(), seq_id, frag_start, frag_end);
		}
	}

	pub fn write_hits(&self, path: std::path::PathBuf) {
		let mut file = File::create(path).expect("Failed to create hits file");
		for (probe, hits) in &self.hits {
			let total_score: f64 = hits.iter().map(|h| h.score).sum();
			if total_score == 0.0 {
				continue;
			}
			for hit in hits {
				writeln!(file, "{}\t{}\t{}\t{}\t{}", probe, hit.seq_id, hit.start, hit.end, hit.score/total_score).unwrap();
			}
		}
	}

}