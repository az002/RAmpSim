use bam::SamReader;
use std::fs::File;
use std::io::BufWriter;
use seq_io::fasta::Reader;

mod utils;
use utils::sim::*;

fn main() {
    let path = "/usr1/aidanz/projects/read_sim/data/alignments/sorted_mock_alignments.sam";
    let mut samreader = SamReader::from_path(path).unwrap();

    let ref_path = "/usr1/aidanz/projects/read_sim/data/fasta/ZymoMOCK/concat_mock.fa";
    let mut ref_reader= Reader::from_path(ref_path).expect("Failed to open reference file");

    let output_path = "/usr1/aidanz/projects/read_sim/data/sim_output/new_reads.fa";
    let output_file = File::create(output_path).expect("Failed to create output file");
    let mut output_writer = BufWriter::new(output_file);

    let read_len = 8000;
    let bucket_size = 1000;

    let mut simulation = Simulation::new(
        read_len,
        bucket_size,
        &mut ref_reader,
        &mut samreader,
        &mut output_writer
    );

    simulation.simulate();
}