use bam::SamReader;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::EncodeUtf16;

fn main() {
    let path = "/usr1/aidanz/projects/read_sim/data/sorted_all_alignments.sam";
    let samreader = SamReader::from_path(path).unwrap();

    let integers_path = "/usr1/aidanz/projects/read_sim/data/seqdump_index.txt";
    let integers_file = File::open(integers_path).expect("Failed to open integers file");
    let reader = BufReader::new(integers_file);

    let mut numbers: Vec<i32> = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if let Ok(num) = line.trim().parse::<i32>() {
            numbers.push(num);
        } else {
            eprintln!("Warning: Could not parse line as integer: {}", line);
        }
    }

    println!("Read {} integers from file", numbers.len());

    for record in samreader {
        let record = record.unwrap();
        let start = record.start();
        let end = record.calculate_end();
        // Find the first integer greater than start using binary search
        let greater_than_start = match numbers.binary_search(&(start as i32)) {
            Ok(exact_idx) => {
                // Found exact match, find the next greater value
                let mut idx = exact_idx;
                while idx < numbers.len() && numbers[idx] <= start as i32 {
                    idx += 1;
                }
                if idx < numbers.len() {
                    Some(numbers[idx])
                } else {
                    None
                }
            },
            Err(insert_idx) => {
                // insert_idx is where the value would be inserted
                if insert_idx < numbers.len() {
                    Some(numbers[insert_idx])
                } else {
                    None
                }
            }
        };

        if let Some(value) = greater_than_start {
            if end >= value{
                println!("skipped: end {} >= value {}", end, value);
                continue;
            }
        }
    }
}