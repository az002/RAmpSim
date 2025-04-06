use bam::SamReader;
use std::fs::File;
use std::io::BufWriter;
use seq_io::fasta::Reader;
use clap::Parser;


mod utils;
use utils::sim::*;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    ref_path: std::path::PathBuf,

    #[arg(short, long)]
    out_path: std::path::PathBuf,

    #[arg(short, long)]
    sam_path: std::path::PathBuf,

    #[arg(short='l', long, default_value_t = 8000)]
    read_len: usize,

    #[arg(short, long, default_value_t = 1000)]
    bucket_size: usize
}

fn main() {
    let args = Args::parse();
    println!("ref_path: {:?}", args.ref_path);
    println!("out_path: {:?}", args.out_path);
    println!("sam_path: {:?}", args.sam_path);
    println!("read_len: {:?}", args.read_len);
    println!("bucket_size: {:?}", args.bucket_size);
    // let path = "/usr1/aidanz/projects/read_sim/data/alignments/sorted_mock_alignments.sam";
    let mut samreader = SamReader::from_path(args.sam_path).unwrap();

    // let ref_path = "/usr1/aidanz/projects/read_sim/data/fasta/ZymoMOCK/concat_mock.fa";
    let mut ref_reader= Reader::from_path(args.ref_path).expect("Failed to open reference file");

    // let output_path = "/usr1/aidanz/projects/read_sim/data/sim_output/new_reads.fa";
    let output_file = File::create(args.out_path).expect("Failed to create output file");
    let mut output_writer = BufWriter::new(output_file);

    // let read_len = 8000;
    // let bucket_size = 1000;

    let mut simulation = Simulation::new(
        args.read_len,
        args.bucket_size,
        &mut ref_reader,
        &mut samreader,
        &mut output_writer
    );

    simulation.simulate();
}