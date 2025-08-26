use core::num;
use std::collections::{VecDeque, HashMap};
use std::hash::Hash;
use bam::{record, Record as BamRecord};
use rand_distr::{Poisson,LogNormal,Distribution};
use rand::Rng;
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


pub trait ProbeTrait {
	fn score(&self) -> f64;
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Probe {
	pub pos_start: usize,
    pub pos_end: usize,
	pub score: f64,
}

impl ProbeTrait for Probe {
	fn score(&self) -> f64 { self.score }
}

pub struct SequenceInput {
	pub samreader: SamReader<BufReader<File>>,
	pub ref_reader: Reader<File>,
	pub output_writer: BufWriter<File>,
}

#[derive(Clone, Default, Debug)]
pub struct Window<T> {
	pub probes : VecDeque<T>,
	center : usize,
	score: f64,
	genome_ind: usize,
}
impl<T: ProbeTrait> Window<T> {

	pub fn get_center(&self) -> Option<&T> { 
		self.probes.get(self.center) 
	}

	pub fn len(&self) -> usize { self.probes.len()}

	pub fn push_back(&mut self, item: T) {
		self.score += item.score();
		self.probes.push_back(item);
	}

	pub fn is_empty(&self) -> bool { self.probes.is_empty() }

	pub fn pop_front(&mut self) -> Option<T> {
		if self.center > 0 {
			self.center -= 1;
		}
		let popped = self.probes.pop_front();
		if popped.is_some() {
			self.score -= popped.as_ref().unwrap().score();
		}
		popped
	}

	pub fn front(&self) -> Option<&T> { self.probes.front() }

	pub fn back(&self) -> Option<&T> { self.probes.back() }

	pub fn set_center(&mut self, index: usize) {
		if index < self.probes.len() {
			self.center = index;
		} else {
			panic!("Index out of bounds for setting center in Window");
		}
	}

	pub fn clear(&mut self) {
		self.probes.clear();
		self.center = 0;
		self.score = 0.0;
	}
}
pub struct Simulator<'a> {
    fragment_len: usize,
	lognorm_sd: f64,
	nfragments: usize,
    use_energy:bool,
	seqinput: &'a mut SequenceInput,
    pub cur_aln: Option<BamRecord>,
    cur_window: Window<Probe>,
    cur_probe: Probe,
	end_reached: bool,
	probe_multiplicities: &'a HashMap<String, usize>,
	vec_windows: Vec<Window<Probe>>,
}

impl<'a> Simulator<'a> {
	pub fn new(fragment_len: usize, lognorm_sd: f64, nfragments: usize, use_energy: bool, seqinput: &'a mut SequenceInput, probe_multiplicities: &'a HashMap<String, usize>) -> Self {
		let mut sim = Simulator {
			fragment_len,
			lognorm_sd,
			nfragments,
			use_energy,
			seqinput,
			cur_aln: None,
			cur_window: Window::default(),
			cur_probe: Probe { pos_start: 0, pos_end: 0, score: 0.0 },
			end_reached: false,
			probe_multiplicities,
			vec_windows: Vec::<Window<Probe>>::new(),
		};
		sim.cur_aln = Simulator::get_next_read(&mut sim.seqinput.samreader);
		if sim.cur_aln.is_some() {
			let aln = sim.cur_aln.as_ref().unwrap();
			sim.cur_probe = Probe {
				pos_start: aln.start() as usize,
				pos_end: (aln.calculate_end()-1) as usize,
				score: Simulator::score_aln(&aln, sim.probe_multiplicities),
			};
		}
		else {
			panic!("No initial alignment found in SAM file");
		}
		sim
	}

	pub fn get_next_read(samreader: &mut SamReader<BufReader<File>>) -> Option<BamRecord> {
        match samreader.next() {
			Some(Ok(aln)) => {
				Some(aln)
			},
			Some(Err(e)) => {
				panic!("Error reading next SAM record: {}", e);
			},
			None => {
				None
			}
        }
    }

	pub fn score_aln(record: &BamRecord, probe_multiplicities: &HashMap<String, usize>) -> f64{
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
		let rname = str::from_utf8(record.name()).unwrap_or("");
		let mult = *probe_multiplicities.get(rname).unwrap_or(&1);
        profile.dot(&SCORE) as f64 * (mult as f64)
	}

	pub fn compute_distr(window: &Window<Probe>, flen: usize) -> Vec<(usize, usize, usize)> {
		let front = window.front().unwrap().pos_start;
		let back = window.back().unwrap().pos_end;

		let start = cmp::min(front, window.get_center().unwrap().pos_end.saturating_sub(flen) + 1);
		let end = cmp::max(back, window.get_center().unwrap().pos_start + flen - 1);

		let mut rectangles: Vec<(usize, usize, usize)> = Vec::new();
		let l = window.len();

		// Slide the window from start to end-flen+1
		let mut window_start = start;
		let mut window_end = window_start + flen - 1;

		let mut probes_in_window: usize = 0;
		let mut probe_starts = Vec::with_capacity(l);
		let mut probe_ends = Vec::with_capacity(l);

		// Collect probe starts and ends
		for probe in &window.probes {
			probe_starts.push(probe.pos_start);
			probe_ends.push(probe.pos_end);
		}

		let mut idx_window_head = 0;
		let mut idx = 0;

		while idx < l
			&& window_start <= probe_starts[idx]
			&& probe_ends[idx] <= window_end
		{
			probes_in_window += 1;
			idx += 1;
		}

		// Store rectangles as (start, width, count)
		while window_end <= end {
			let s_diff = if probes_in_window > 0 { probe_starts[idx_window_head].saturating_sub(window_start) } else { usize::MAX };
			let e_diff = if idx < l { probe_ends[idx].saturating_sub(window_end) } else { usize::MAX };
			let max_jump = end.saturating_sub(window_end);
        	let jump = cmp::min(cmp::min(s_diff, e_diff), max_jump);

			rectangles.push((window_start, jump+1, cmp::max(probes_in_window, 1)));
			window_start += jump + 1;
			window_end = window_start + flen - 1;

			// Remove probes that are no longer in the window
			while idx_window_head < l && probe_ends[idx_window_head] < window_start {
				probes_in_window -= 1;
				idx_window_head += 1;
			}

			while idx < l
				&& window_start <= probe_starts[idx]
				&& probe_ends[idx] <= window_end
			{
				probes_in_window += 1;
				idx += 1;
			}
		}

		rectangles
	}

	pub fn simulate_telseq_window_fragments(
        window: &Window<Probe>, 
        count: usize, 
        score: f64, 
        flen: usize,
		lognorm_sd: f64,
        nfragments: usize,
        use_energy: bool,
        total_score: f64, 
        total_count: usize, 
        seq: &[u8], 
        cur_genome: &[u8], 
        output_writer: &mut BufWriter<File>,
		split: f64,
    ) -> usize{
        let mut lambda = count as f64/ total_count as f64 * split * nfragments as f64;
        if use_energy {
            lambda = score / total_score * split * nfragments as f64;
        }
		let rectangles = Simulator::compute_distr(window, flen);
		let poisson = Poisson::new(lambda).unwrap();
		let lognorm = LogNormal::new((flen as f64).ln(), lognorm_sd).unwrap();
		let n: usize = poisson.sample(&mut rand::rng()) as usize;

		let mut rng = rand::rng();

		// If there are no valid rectangles, return early
		if rectangles.is_empty() {
			return 0;
		}

		// Compute the total weight and cumulative weights for O(1) sampling
		let mut cum_weights = Vec::with_capacity(rectangles.len());
		let mut total_weight = 0.0;
		for 	&(_, width, count) in &rectangles {
			total_weight += (width * count) as f64;
			cum_weights.push(total_weight);
		}

		// Sample fragment start positions
		for _ in 0..n{
			let x = rng.random_range(0.0..total_weight);
			// Binary search to find the rectangle
			let rect_idx = cum_weights.binary_search_by(|w| w.partial_cmp(&x).unwrap()).unwrap_or_else(|i| i);
			let &(rect_start, rect_width, _) = &rectangles[rect_idx];
			// Uniformly pick a position within the rectangle
			let offset = rng.random_range(0..rect_width);
			let frag_start = rect_start + offset;
			let frag_end = (frag_start + lognorm.sample(&mut rng) as usize).min(seq.len());
			// println!("{} {} {} {}", frag_start, frag_end, str::from_utf8(cur_genome).unwrap_or("genome"), seq.len()); // Debugging output
			let frag_seq = &seq[frag_start..frag_end];
			writeln!(
				output_writer,
				">{}_{}_{}\n{}",
				std::str::from_utf8(cur_genome).unwrap_or("genome"),
				frag_start,
				frag_end,
				std::str::from_utf8(frag_seq).unwrap_or("")
			).unwrap();
			output_writer.flush().unwrap();
		}
		n
	}

	pub fn compute_windows(&mut self) -> (f64, usize) {
		let mut genome_ind = -1;
		let mut total_score = 0.0;
        let mut total_count: usize = 0;
		loop {
			if self.end_reached {
				break;
			}
			genome_ind += 1;
			self.cur_window.genome_ind = genome_ind as usize;

			loop {

				if self.cur_aln.is_none() {
					self.end_reached = true;
					break;
				}

				let aln = self.cur_aln.as_ref().unwrap();

				if aln.ref_id() != genome_ind {
					if !self.cur_window.is_empty() {
						self.vec_windows.push(self.cur_window.clone());
						total_score += self.cur_window.score;
                        total_count += self.cur_window.len();
						self.cur_window.clear();
						// Generate fragments for the current window before moving to the next genome
					}
					break; // Skip to the next genome if the current reference ID does not match
				}

				if self.cur_window.is_empty() || self.cur_probe.pos_end <= self.cur_window.get_center().unwrap().pos_start + self.fragment_len - 1 {
					self.cur_window.push_back(self.cur_probe);
				} else {
					self.vec_windows.push(self.cur_window.clone());
					total_score += self.cur_window.score;
                    total_count += self.cur_window.len();
					self.cur_window.push_back(self.cur_probe);
					self.cur_window.set_center(self.cur_window.len() - 1);

					while !self.cur_window.is_empty() && self.cur_window.front().unwrap().pos_start < self.cur_probe.pos_end - self.fragment_len + 1 {
						self.cur_window.pop_front();
					}
				}

				self.cur_aln = Simulator::get_next_read(&mut self.seqinput.samreader);
					if self.cur_aln.is_some() {
						let aln = self.cur_aln.as_ref().unwrap();
						self.cur_probe = Probe {
							pos_start: aln.start() as usize,
							pos_end: (aln.calculate_end()-1) as usize,
							score: Simulator::score_aln(aln, self.probe_multiplicities),
						};
					} 				
			}
		}
		(total_score, total_count)
	}
	pub fn simulate_telseq_fragments(&mut self, split: f64) -> usize {
		let (total_score, total_count) = self.compute_windows();
		let mut cur_genome_ind = 0;
		let mut ind = 0;
		let mut fragments_generated = 0;
		while let Some(record) = self.seqinput.ref_reader.next() {
			let record = record.expect("Error reading record");
			let cur_genome = record.id_bytes();
			let seq = record.seq();
			
			while ind < self.vec_windows.len() && self.vec_windows[ind].genome_ind == cur_genome_ind {
				let window = &self.vec_windows[ind];
				fragments_generated += Simulator::simulate_telseq_window_fragments(window, window.len(), window.score, self.fragment_len, self.lognorm_sd, self.nfragments, self.use_energy, total_score, total_count, seq, cur_genome, &mut self.seqinput.output_writer, split);
				ind += 1;
			}
			cur_genome_ind += 1;
		}
		fragments_generated
	}

	pub fn simulate_genome_background(
		seqs: &Vec<Vec<u8>>,
		cur_genome: String,
		flen: usize,
		lognorm_sd: f64,
		nfrags: usize,
		output_writer: &mut BufWriter<File>,
	) {
		let lognorm = LogNormal::new((flen as f64).ln(), lognorm_sd).unwrap();
		let total_length = seqs.iter().map(|seq| seq.len()).sum::<usize>();
		let mut rng = rand::rng();

		for seq in seqs {
			let n = seq.len()/total_length * nfrags;
			for _ in 0..n {
				let frag_start = rng.random_range(0..seq.len());
				let frag_end = frag_start + lognorm.sample(&mut rng) as usize;
				if frag_end <= seq.len() {
					let frag_seq = &seq[frag_start..frag_end];
					writeln!(
						output_writer,
						">{}_{}_{}\n{}",
						cur_genome,
						frag_start,
						frag_end,
						std::str::from_utf8(frag_seq).unwrap_or("")
					).unwrap();
				}
			}
		}
	}
	

	pub fn simulate_background_fragments(
		&mut self,
		nfrags: usize,
		abundances: &HashMap<String, f64>,
		seqid2refid: &HashMap<String, String>,
	) {
		let nfrags = self.nfragments - nfrags;
		let mut seqs = Vec::new();
		let record = self.seqinput.ref_reader.next().expect("Error reading first record").unwrap();
		seqs.push(record.seq().to_vec());
		let mut cur_genome = seqid2refid.get(str::from_utf8(record.id_bytes()).unwrap()).expect("Failed to get sequence ID").to_string();
		while let Some(record) = self.seqinput.ref_reader.next() {
			let record = record.expect("Error reading record");
			let id = seqid2refid.get(str::from_utf8(record.id_bytes()).unwrap()).expect("Failed to get sequence ID").to_string();
			if *id == *cur_genome {
				seqs.push(record.seq().to_vec());
			}
			else {
				let n = abundances.get(&id).unwrap() * nfrags as f64;
				Simulator::simulate_genome_background(&seqs, cur_genome, self.fragment_len, self.lognorm_sd, n as usize, &mut self.seqinput.output_writer);
				seqs.clear();
				cur_genome = id;
				seqs.push(record.seq().to_vec());
			}
		}
		let n = abundances.get(&cur_genome).unwrap() * nfrags as f64;
		Simulator::simulate_genome_background(&seqs, cur_genome, self.fragment_len, self.lognorm_sd, n as usize, &mut self.seqinput.output_writer);
	}

	pub fn simulate(
		&mut self,
		abundances: &HashMap<String, f64>,
		seqid2refid: &HashMap<String, String>,
		split: f64,
	) {
		self.seqinput.ref_reader.next();
		let pos1 = self.seqinput.ref_reader.position().unwrap().clone();
		let _ = self.seqinput.ref_reader.seek(&pos1);
		let fragments_generated = self.simulate_telseq_fragments(split);
		let _ = self.seqinput.ref_reader.seek(&pos1);
		self.simulate_background_fragments(fragments_generated, abundances, seqid2refid);
	}
}