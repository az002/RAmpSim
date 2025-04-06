use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::io::BufReader;
use bam::SamReader;
use bam::Record as BamRecord;
use seq_io::fasta::{Reader, Record};
use std::cmp::min;
use rand_distr::{Poisson,Distribution};
use rand::Rng;

pub struct Simulator {
    read_len : usize,
    bucket_size : usize,
    pub aln : BamRecord,
    cur_aln_start : usize,
    cur_bucket_start : usize,
    cur_bucket_end : usize,
    cur_count : usize,
    end_reached : bool,
    nreads : usize
}

impl Simulator {
    pub fn new(read_len: usize, bucket_size: usize) -> Simulator {
        Simulator {
            read_len,
            bucket_size,
            aln: BamRecord::new(),
            cur_aln_start: 0,
            cur_bucket_start: 0,
            cur_bucket_end: 0,
            cur_count: 0,
            end_reached: false,
            nreads: 0
        }
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

    pub fn generate_reads(
        &mut self,
        output_writer: &mut BufWriter<File>,
        read_len: usize,
        score: usize,
        seq: &[u8],
        cur_bucket_start: usize,
        cur_bucket_end: usize,
        nreads: usize,
        cur_genome: &[u8],
    ) -> usize {
        let poisson = Poisson::new(score as f64).unwrap();
        let s = poisson.sample(&mut rand::rng());
        println!(
            "Generating {} reads for bucket {} to {} containing {} reads in genome {}",
            s, cur_bucket_start, cur_bucket_end, score, str::from_utf8(cur_genome).unwrap()
        );
        let bucket_size = cur_bucket_end - cur_bucket_start + 1;
        let valid_start = cur_bucket_start.saturating_sub(read_len-bucket_size);
        let valid_end = cur_bucket_start;

        // Generate reads if we have a valid range
        if valid_start < valid_end {
            let mut rng = rand::rng();
            for i in 0..s as usize {
                // Generate a random starting position
                let start_pos = rng.random_range(valid_start..=valid_end);
                let end_pos = min(start_pos + read_len, seq.len());
                
                // Extract the read from the sequence
                let read = &seq[start_pos..end_pos];
                // Write the read to the output file
                let id = format!("read_{}_{}_{}_{}", nreads+i, str::from_utf8(cur_genome).unwrap(), start_pos, end_pos);
                let _ = seq_io::fasta::write_parts(&mut *output_writer, id.as_bytes(), None, read);
                
                // Do something with the read (save to file, store in collection, etc.)
                //println!("Generated read at position {}", start_pos+cur_start);
            }
        }
        nreads + (s as usize)
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
        read_len: usize,
        bucket_size: usize,
        ref_reader: &'a mut Reader<File>,
        samreader: &'a mut SamReader<BufReader<File>>,
        output_writer: &'a mut BufWriter<File>,
    ) -> Simulation<'a> {
        Simulation {
            sim: Simulator::new(
                read_len,
                bucket_size,
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
            }
            let mut cur_bucket_end = min(self.sim.cur_bucket_start + self.sim.bucket_size - 1, seq.len()-1);
            loop {
                if genome_ind == self.sim.aln.ref_id() && self.sim.cur_aln_start <= cur_bucket_end {
                    self.sim.cur_count += 1;
                    if self.sim.get_next_read(self.samreader) {
                        self.sim.end_reached = true;
                        let _ = self.sim.generate_reads(self.output_writer, self.sim.read_len, self.sim.cur_count, seq, self.sim.cur_bucket_start, cur_bucket_end, self.sim.nreads, cur_genome);
                        break;
                    }

                    self.sim.cur_aln_start = self.sim.aln.start() as usize;
                    continue;
                }

                let _ = self.sim.generate_reads(self.output_writer, self.sim.read_len, self.sim.cur_count, seq, self.sim.cur_bucket_start, cur_bucket_end, self.sim.nreads, cur_genome);
                self.sim.cur_bucket_start = self.sim.cur_aln_start;
                cur_bucket_end = min(self.sim.cur_bucket_start + self.sim.bucket_size - 1, seq.len()-1);
                self.sim.cur_count = 0;
                if self.sim.aln.ref_id() != genome_ind { break; }
            }
        }
    }
}

