use bam::SamReader;
use std::fs::File;
use std::io::{BufWriter, BufRead};
use std::collections::HashMap;
use seq_io::fasta::{Reader, Record};
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

    #[arg(short='l', long="flen", default_value_t = 8000)]
    fragment_len: usize,

    #[arg(short,long)]
    temperature: f64,

    #[arg(short, long="mult-file")]
    probe_multiplicity: Option<std::path::PathBuf>,

    #[arg(long="seqids")]
    sequence_ids: std::path::PathBuf,

    #[arg(short,long)]
    abundances: std::path::PathBuf,

    #[arg(short,long="nfrag")]
    nfragments: usize,

    #[arg(long)]
    split: f64,

    #[arg(long, default_value_t = 0.4)]
    lognorm_sd: f64
}

fn main() {
    let args = Args::parse();
    eprintln!("ref_path: {:?}", args.ref_path);
    eprintln!("out_path: {:?}", args.out_path);
    eprintln!("sam_path: {:?}", args.sam_path);
    eprintln!("fragment_len: {:?}", args.fragment_len);
    eprintln!("nfragments: {:?}", args.nfragments);
    eprintln!("temperature: {:?}", args.temperature);

    let mut probe_alignments = SamReader::from_path(args.sam_path).unwrap();
    let mut output = BufWriter::new(File::create(args.out_path).expect("Failed to create output file"));

    let probe_multiplicities: HashMap<String, usize> = match args.probe_multiplicity {
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

    let abundances: HashMap<String, f64> = {
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

    let seqid2refid: HashMap<String, String> = {
        let file = File::open(args.sequence_ids).expect("Failed to open sequence IDs file");
        let reader = std::io::BufReader::new(file);
        let mut map: HashMap<String, String> = HashMap::new();
        for line in reader.lines() {
            let line = line.expect("Failed to read line");
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 2 {
                let id = parts[0].to_string();
                let seq_id = parts[1].to_string();
                map.insert(id, seq_id);
            }
        }
        map
    };

    let mut ref_reader = Reader::from_path(args.ref_path).expect("Failed to open reference file");
    let ref_seqs: HashMap<String, String> = ref_reader
        .records()
        .map(|r| {
            let record = r.expect("Failed to read record");
            let id = record.id().unwrap().to_string();
            let seq = String::from_utf8(record.seq().to_vec()).expect("Failed to convert sequence to string");
            (id, seq)
        })
        .collect();


    // println!("Probe multiplicities: {:?}", probe_multiplicities);
    // println!("Abundances: {:?}", abundances);
    // println!("Sequence ID to Reference ID mapping: {:?}", sequenceid2refid);
    // let read_len = 8000;
    // let bucket_size = 1000;

    // simulation.sim.score_aln();

    let mut sim = Simulator::new(
        args.fragment_len, 
        args.lognorm_sd, 
        args.nfragments,
        args.split,
        args.temperature,
        probe_multiplicities,
        abundances,
        seqid2refid,
    );
    sim.compute_telseq_distr(&mut probe_alignments);
    sim.sample_telseq(&mut output, &ref_seqs);
    sim.sample_background(&mut output, &ref_seqs);
    sim.write_hits("/usr1/aidanz/projects/read_sim/data/sim_output/hits.txt");
	// sim.simulate(args.split);
}