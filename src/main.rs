use bam::SamReader;
use bam::Record as BamRecord;
#[allow(unused_imports)]
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use seq_io::fasta::{Reader, Record};
use std::cmp::min;
use rand_distr::{Poisson,Distribution};
use rand::Rng;

#[allow(unused)]

fn generate_reads(mut writer: &mut BufWriter<File>, read_len: usize, score: i32, seq: &[u8], bucket_start:usize, bucket_end:usize, nreads:usize, genome:&[u8]) -> usize {
    let poisson = Poisson::new(score as f64).unwrap();
    let s = poisson.sample(&mut rand::rng());
    // Calculate valid start position range
    let bucket_size = bucket_end - bucket_start + 1;
    let valid_start = bucket_start.saturating_sub(read_len-bucket_size);
    let valid_end = bucket_start;

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
            let id = format!("read_{}_{}_{}_{}", nreads+i, str::from_utf8(genome).unwrap(), start_pos, end_pos);
            seq_io::fasta::write_parts(&mut writer, id.as_bytes(), None, read);
            
            // Do something with the read (save to file, store in collection, etc.)
            //println!("Generated read at position {}", start_pos+cur_start);
        }
    }
    nreads + (s as usize)
}

fn read_next(samreader: &mut SamReader<BufReader<File>>, align : &mut BamRecord) -> bool {
    match samreader.next() {
        Some(Ok(aln)) => {
            *align = aln;
            false
        },
        _ => {
            true
        }
    }
}

fn main() {
    let path = "/usr1/aidanz/projects/read_sim/data/alignments/sorted_mock_alignments.sam";
    let mut samreader = SamReader::from_path(path).unwrap();

    let ref_path = "/usr1/aidanz/projects/read_sim/data/fasta/ZymoMOCK/concat_mock.fa";
    let mut ref_reader= Reader::from_path(ref_path).expect("Failed to open reference file");

    let output_path = "/usr1/aidanz/projects/read_sim/data/sim_output/reads.fa";
    let output_file = File::create(output_path).expect("Failed to create output file");
    let mut output_writer = BufWriter::new(output_file);

    let read_len = 8000;
    let bucket_size = 1000;

    let mut cur_aln = BamRecord::new();
    read_next(&mut samreader, &mut cur_aln);

    let mut cur_aln_start = cur_aln.start() as usize;
    //let mut cur_aln_end = cur_aln.calculate_end() as usize;

    let mut cur_bucket_start = cur_aln.start() as usize;
    // let mut cur_start = 0;
    let mut genome_ind = -1;
    let mut cur_count = 0;
    let mut end_reached = false;
    let mut nreads = 0;

    for record in ref_reader.records() {
        if end_reached { break; }

        let record = record.expect("Error reading record");
        let cur_genome = record.id_bytes();
        let seq = record.seq();
        genome_ind += 1;

        if cur_aln.ref_id() != genome_ind {
            continue;
        }
        let mut cur_bucket_end = min(cur_bucket_start + bucket_size - 1, seq.len()-1);
        loop {
            if genome_ind == cur_aln.ref_id() && cur_aln_start <= cur_bucket_end {
                cur_count += 1;
                if read_next(&mut samreader, &mut cur_aln) {
                    end_reached = true;
                    break;
                }

                cur_aln_start = cur_aln.start() as usize;
                continue;
            }

            nreads = generate_reads(&mut output_writer, read_len, cur_count, seq, cur_bucket_start, cur_bucket_end, nreads, cur_genome);
            println!("Generated {} reads for genome {} bucket start: {}, bucket end: {}", cur_count, genome_ind, cur_bucket_start, cur_bucket_end);
            cur_bucket_start = cur_aln_start;
            cur_bucket_end = min(cur_bucket_start + bucket_size - 1, seq.len()-1);
            cur_count = 0;
            if cur_aln.ref_id() != genome_ind { break; }
        }
    }
}