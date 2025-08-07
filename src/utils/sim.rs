use std::fs::File;
use std::io::BufWriter;
use std::io::BufReader;
use bam::SamReader;
use bam::Record as BamRecord;
use nalgebra::ArrayStorage;
use nalgebra::SVector;
use seq_io::fasta::{Reader, Record};
use std::cmp::min;
use rand_distr::{Poisson,Distribution};
use rand::Rng;
use std::collections::HashMap;
use lazy_static::lazy_static;

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

lazy_static! {
    static ref TRIMER_INDEX: HashMap<u8, usize> = HashMap::from([
        (0, 0), (63, 0), (1, 1), (47, 1),
        (2, 2), (31, 2), (3, 3), (15, 3),
        (4, 4), (59, 4), (5, 5), (43, 5),
        (6, 6), (27, 6), (7, 7), (11, 7),
        (8, 8), (55, 8), (9, 9), (39, 9),
        (10, 10), (23, 10), (12, 11), (51, 11),
        (13, 12), (35, 12), (14, 13), (19, 13),
        (16, 14), (62, 14), (17, 15), (46, 15),
        (18, 16), (30, 16), (20, 17), (58, 17),
        (21, 18), (42, 18), (22, 19), (26, 19),
        (24, 20), (54, 20), (25, 21), (38, 21),
        (28, 22), (50, 22), (29, 23), (34, 23),
        (32, 24), (61, 24), (33, 25), (45, 25),
        (36, 26), (57, 26), (37, 27), (41, 27),
        (40, 28), (53, 28), (44, 29), (49, 29),
        (48, 30), (60, 30), (52, 31), (56, 31),
    ]);
}

const SCORE : SVector<f64, 32> = SVector::<f64,32>::from_array_storage(ArrayStorage([[-1.0333523772924682,  -1.458111012202453,  -1.035949547554549,  -0.7124736682600892, 
    -1.8355462789826653,  -3.2816907774428596,  -2.9387615700865912,  -2.4322581762379056, 
    -6.050351136676882e-26,  -1.16892836679913,  -1.2733818523786504,  -0.42044691866600825, 
    -2.2817669581121254,  -1.4822900733089985,  -0.9125362005565887,  -3.670822543027401e-23, 
    -1.273598247757921,  -1.3530081563908097e-28,  -0.19137969659843326,  -1.2926131142356785, 
    -0.3257778299362714,  -2.1797369472925356,  -1.3319137918104318,  -1.7767543614895418, 
    -2.617678120625974,  -1.5216318279285466e-24,  -0.8229383115364808,  -2.1347399492841697, 
    -1.502770970281663e-38,  -0.8507439536421594,  -0.31680823162220356,  -0.06611591220578053]]));

pub struct Simulator {
    fragment_len : usize,
    bucket_size : usize,
    use_energy: bool,
    pub aln : BamRecord,
    cur_aln_start : usize,
    cur_bucket_start : usize,
    cur_bucket_end : usize,
    cur_count : usize,
    cur_score: f64,
    end_reached : bool,
    nreads : usize
}

impl Simulator {
    pub fn new(fragment_len: usize, bucket_size: usize, use_energy: bool) -> Simulator {
        Simulator {
            fragment_len,
            bucket_size,
            use_energy,
            aln: BamRecord::new(),
            cur_aln_start: 0,
            cur_bucket_start: 0,
            cur_bucket_end: 0,
            cur_count: 0,
            cur_score: 0.,
            end_reached: false,
            nreads: 0
        }
    }

    pub fn score_aln(&mut self) {
        let mut trimer = 0u8;
        let mask:u8 = 0b00111111;
        let mut profile = SVector::<f64, 32>::zeros();
        let mut count : u8 = 0;
        for entry in self.aln.alignment_entries().unwrap() {
            if entry.is_seq_match() {
                unsafe {
                    trimer = (trimer << 2) & mask | (*BYTE2BITS.as_ptr().add(entry.record_nt().unwrap() as usize));
                }
                count += 1;
                if count >= 3 {
                    let mut index = 0;
                    if let Some(&i) = TRIMER_INDEX.get(&trimer) {
                        index = i;
                    }
                    profile[index] += 1.;
                }
            }
            else {
                trimer = 0;
                count = 0;
            }
        }
        self.cur_score -= profile.dot(&SCORE) as f64;
    }

    pub fn init(&mut self, samreader: &mut SamReader<BufReader<File>>) {
        self.get_next_read(samreader);
        self.cur_aln_start = self.aln.start() as usize;
        self.cur_bucket_start = self.aln.start() as usize;
    }

    pub fn get_next_read(&mut self, samreader: &mut SamReader<BufReader<File>>) -> bool {
        match samreader.next() {
            Some(Ok(aln)) => {
                self.aln = aln;
                false
            },
            _ => {
                true
            }
        }
    }

    pub fn generate_fragments(
        &mut self,
        output_writer: &mut BufWriter<File>,
        seq: &[u8],
        cur_genome: &[u8],
    ) {
        let lambda: f64;
        if self.use_energy {
            lambda = self.cur_score/20.;
        } else {
            lambda = self.cur_count as f64;
        }
        let poisson = Poisson::new(lambda).unwrap();
        let s = poisson.sample(&mut rand::rng());
        println!(
            "Generating {} fragments for bucket [{} - {}] containing {} probes with score {} in genome {}",
            s, self.cur_bucket_start, self.cur_bucket_end, self.cur_count, self.cur_score, str::from_utf8(cur_genome).unwrap()
        );
        let bucket_size = self.cur_bucket_end - self.cur_bucket_start + 1;
        let valid_start = self.cur_bucket_start.saturating_sub(self.fragment_len-bucket_size);
        let valid_end = self.cur_bucket_start;

        // Generate reads if we have a valid range
        if valid_start < valid_end {
            let mut rng = rand::rng();
            for i in 0..s as usize {
                // Generate a random starting position
                let start_pos = rng.random_range(valid_start..=valid_end);
                let end_pos = min(start_pos + self.fragment_len, seq.len());
                
                // Extract the read from the sequence
                let read = &seq[start_pos..end_pos];
                // Write the read to the output file
                let id = format!("read_{}_{}_{}_{}", self.nreads+i, str::from_utf8(cur_genome).unwrap(), start_pos, end_pos);
                let _ = seq_io::fasta::write_parts(&mut *output_writer, id.as_bytes(), None, read);
                
                // Do something with the read (save to file, store in collection, etc.)
                //println!("Generated read at position {}", start_pos+cur_start);
            }
        }
        self.nreads += s as usize;
    }

}

pub struct Simulation<'a> {
    pub sim : Simulator,
    ref_reader : &'a mut Reader<File>,
    pub samreader : &'a mut SamReader<BufReader<File>>,
    output_writer : &'a mut BufWriter<File>,
}

impl<'a> Simulation<'a> {
    pub fn new(
        fragment_len: usize,
        bucket_size: usize,
        use_energy: bool,
        ref_reader: &'a mut Reader<File>,
        samreader: &'a mut SamReader<BufReader<File>>,
        output_writer: &'a mut BufWriter<File>,
    ) -> Simulation<'a> {
        Simulation {
            sim: Simulator::new(
                fragment_len,
                bucket_size,
                use_energy
            ),
            ref_reader,
            samreader,
            output_writer,
        }
    }

    pub fn simulate(&mut self) {
        let mut genome_ind = -1;
        self.sim.init(self.samreader);

        for record in self.ref_reader.records() {
            if self.sim.end_reached { break; }

            let record = record.expect("Error reading record");
            let cur_genome = record.id_bytes();
            let seq = record.seq();
            genome_ind += 1;

            if self.sim.aln.ref_id() != genome_ind {
                continue;
                // Skip to the next genome if the current reference ID does not match
            }

            self.sim.cur_bucket_end = min(self.sim.cur_bucket_start + self.sim.bucket_size - 1, seq.len()-1);
            loop {
                if genome_ind == self.sim.aln.ref_id() && self.sim.cur_aln_start <= self.sim.cur_bucket_end {
                    self.sim.cur_count += 1;
                    self.sim.score_aln();
                    if self.sim.get_next_read(self.samreader) {
                        self.sim.end_reached = true;
                        self.sim.generate_fragments(self.output_writer, seq, cur_genome);
                        break;
                    }

                    self.sim.cur_aln_start = self.sim.aln.start() as usize;
                    continue;
                }

                self.sim.generate_fragments(self.output_writer, seq, cur_genome);
                self.sim.cur_bucket_start = self.sim.cur_aln_start;
                self.sim.cur_bucket_end = min(self.sim.cur_bucket_start + self.sim.bucket_size - 1, seq.len()-1);
                self.sim.cur_count = 0;
                self.sim.cur_score = 0.;
                if self.sim.aln.ref_id() != genome_ind { break; }
            }
        }
    }
}

