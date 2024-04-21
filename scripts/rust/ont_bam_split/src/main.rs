extern crate bam;
extern crate clap;
extern crate levenshtein;

use std::fs::File;
use std::io::{BufWriter, stdout, Write};
use bam::{Reader, Record};
use clap::{Arg, App};
use levenshtein::levenshtein;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = App::new("My Program")
        .version("1.0")
        .author("Your Name <your.email@example.com>")
        .about("Does awesome things")
        .arg(Arg::with_name("INPUT")
            .help("Sets the input file to use")
            .required(true)
            .index(1))
        .arg(Arg::with_name("ADAPTER")
            .short("a")
            .long("adapter")
            .value_name("ADAPTER_SEQUENCE")
            .help("Sets the adapter sequence to use")
            .takes_value(true))
        .arg(Arg::with_name("ERROR")
            .short("e")
            .long("error")
            .value_name("ERROR_RATE")
            .help("Sets the maximum allowed error rate for adapter sequence matching")
            .takes_value(true))
        .arg(Arg::with_name("OUTPUT")
            .short("o")
            .long("output")
            .value_name("OUTPUT_PATH")
            .help("Sets the output file path")
            .takes_value(true))
        .get_matches();

    let input_file = File::open(matches.value_of("INPUT").unwrap())?;
    let mut reader = Reader::new(input_file);

    let output: Box<dyn Write> = match matches.value_of("OUTPUT") {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(stdout()),
    };
    let mut writer = bam::Writer::from_stream(BufWriter::new(output), &reader.header().clone())?;

    let adapter_sequence = matches.value_of("ADAPTER").unwrap();
    let max_error_rate = matches.value_of("ERROR").unwrap().parse::<usize>()?;

    for result in reader.records() {
        let mut record = result?;

        let seq = record.seq().to_owned();

        // Calculate the Levenshtein distance between the sequence and the adapter sequence
        let distance = levenshtein(seq, adapter_sequence);

        // Check if the distance is less than or equal to the maximum allowed error rate
        if distance <= max_error_rate {
            let pos = seq.find(adapter_sequence).unwrap();
            let (left, right) = seq.split_at(pos);
            record.set_seq(right)?;
            record.tags_mut().push_str(b"LT", left)?;
        }

        writer.write(&record)?;
    }

    Ok(())
}
