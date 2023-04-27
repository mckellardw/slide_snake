//// command line to build and run:
// cd scripts/rust/h5_to_list_rust
// cargo build --release
// cd ..
// ./h5_to_list_rust/target/release/h5_to_list_rust input_file.h5 output_file.txt

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;

use hdf5::{H5Type, Dataset, Group, Reader};
use itertools::izip;

fn seq_decode(seq_int: u64, seq_len: usize) -> String {
    const ATCG_BASES: [char; 4] = ['A', 'C', 'T', 'G'];  // Changed the order of bases
    let mut seqs = String::with_capacity(seq_len);
    for i in 0..seq_len {
        let tint = (seq_int >> (i * 2)) & 3;
        seqs.push(ATCG_BASES[tint as usize]);
    }
    seqs
}

fn process_group(group: &Group, output_file: Arc<Mutex<File>>, seq_length: usize, n_threads: usize, verbose: bool) {
    for key in group.member_names().unwrap() {
        let value = group.dataset(&key).unwrap();
        let dataset_array: Vec<u64> = value.read_1d().unwrap();

        let chunk_size = dataset_array.len() / n_threads;
        let mut thread_handles = vec![];

        for chunk_start in (0..dataset_array.len()).step_by(chunk_size) {
            let chunk_end = (chunk_start + chunk_size).min(dataset_array.len());
            let chunk = dataset_array[chunk_start..chunk_end].to_vec();
            let output_file_clone = Arc::clone(&output_file);

            let handle = thread::spawn(move || {
                let mut local_results = Vec::new();
                for (row_idx, row) in chunk.iter().enumerate() {
                    let nucleotide_pairs = seq_decode(*row, seq_length);
                    let result = format!("{}\t{}\t{}\n", nucleotide_pairs, row_idx, chunk_start);
                    local_results.push(result);
                }
                let mut file = output_file_clone.lock().unwrap();
                for result in local_results {
                    write!(file, "{}", result).unwrap();
                }
            });
            thread_handles.push(handle);
        }

        for handle in thread_handles {
            handle.join().unwrap();
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: h5_to_list input_file output_file");
        std::process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];
    let seq_length = 25;
    let n_threads = 4;
    let verbose = false;

    let file = hdf5::File::open(input_file).unwrap();
    let out_file = Arc::new(Mutex::new(File::create(output_file).unwrap()));
    let group = file.group("/").unwrap();

    process_group(&group, out_file, seq_length, n_threads, verbose);
}
