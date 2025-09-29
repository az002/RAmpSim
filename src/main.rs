use bam::SamReader;
use std::fs::File;
use std::io::{BufWriter, BufRead};
use std::collections::HashMap;
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
    fragment_len: usize,

    #[arg(short, long, action)]
    use_energy: bool,

    #[arg(short, long)]
    probe_multiplicity: Option<std::path::PathBuf>,

    #[arg(long)]
    sequence_ids: std::path::PathBuf,

    #[arg(short,long)]
    abundances: std::path::PathBuf,

    #[arg(short,long)]
    nfragments: usize,

    #[arg(long)]
    pass_num: usize,

    #[arg(long)]
    split: f64,

    #[arg(long, default_value_t = 0.4)]
    lognorm_sd: f64
}

fn main() {
    let args = Args::parse();
    println!("ref_path: {:?}", args.ref_path);
    println!("out_path: {:?}", args.out_path);
    println!("sam_path: {:?}", args.sam_path);
    println!("fragment_len: {:?}", args.fragment_len);
    println!("nfragments: {:?}", args.nfragments);
    println!("pass_num: {:?}", args.pass_num);
    println!("use_energy: {:?}", args.use_energy);
    // let path = "/usr1/aidanz/projects/read_sim/data/alignments/sorted_mock_alignments.sam";

    // let ref_path = "/usr1/aidanz/projects/read_sim/data/fasta/ZymoMOCK/concat_mock.fa";

    // let output_path = "/usr1/aidanz/projects/read_sim/data/sim_output/new_reads.fa";

    //loop through sam file once to get number of times each probe was aligned
    let records = SamReader::from_path(&args.sam_path).unwrap();
    let mut probe_counts: HashMap<String, usize> = HashMap::new();
    for record in records {
        let record = record.unwrap();
        let probe = String::from_utf8(record.name().to_vec()).unwrap();
        *probe_counts.entry(probe).or_insert(0) += 1;
    }

    let mut seq_in = SequenceInput {
		samreader: SamReader::from_path(args.sam_path).unwrap(),
		ref_reader: Reader::from_path(args.ref_path).expect("Failed to open reference file"),
		output_writer: BufWriter::new(File::create(args.out_path).expect("Failed to create output file")),
	};

    let probe_multiplicities = match args.probe_multiplicity {
        Some(path) => {
            let file = File::open(path).expect("Failed to open probe multiplicity file");
            let reader = std::io::BufReader::new(file);
            let mut map = std::collections::HashMap::new();
            for line in reader.lines() {
                let line = line.expect("Failed to read line");
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() == 2 {
                    let probe = parts[0].to_string();
                    let count: usize = parts[1].parse().expect("Failed to parse count");
                    map.insert(probe, count);
                }
            }
            map
        },
        None => std::collections::HashMap::new(),
    };

    let abundances = {
        let file = File::open(args.abundances).expect("Failed to open abundances file");
        let reader = std::io::BufReader::new(file);
        let mut map = std::collections::HashMap::new();
        for line in reader.lines() {
            let line = line.expect("Failed to read line");
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 2 {
                let id = parts[0].to_string();
                let count: f64 = parts[1].parse().expect("Failed to parse count");
                map.insert(id, count);
            }
        }
        map
    };
    
    let sequenceid2refid = {
        let file = File::open(args.sequence_ids).expect("Failed to open sequence IDs file");
        let reader = std::io::BufReader::new(file);
        let mut map = HashMap::new();
        for line in reader.lines() {
            let line = line.expect("Failed to read line");
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 2 {
                let id = parts[0].to_string();
                let seq_id: String = parts[1].to_string();
                map.insert(id, seq_id);
            }
        }
        map
    };

    // println!("Probe multiplicities: {:?}", probe_multiplicities);
    // println!("Abundances: {:?}", abundances);
    // println!("Sequence ID to Reference ID mapping: {:?}", sequenceid2refid);
    // let read_len = 8000;
    // let bucket_size = 1000;

    // simulation.sim.score_aln();

    let mut sim = Simulator::new(args.fragment_len, args.lognorm_sd, args.nfragments/args.pass_num, args.use_energy, &mut seq_in, &probe_multiplicities, probe_counts, abundances, sequenceid2refid);
	sim.simulate(args.split);
}