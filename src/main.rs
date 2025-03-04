#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
mod bit_tree;
mod lcsk;

use poa::*;
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::{thread, fs::File, io::{self, BufRead}, path::Path};
use lcsk::{lcsk_pipeline, threaded_lcsk_pipeline};
use std::time::Instant;

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -2;
const GAP_OPEN_SCORE: i32 = 2;

fn main() {
    bench_mark_all ();
}

fn bench_mark_all () {
    let test_cases = vec![
        ("./data/Synthetic10.fa", 5, 20, 2),
        //("./data/Synthetic100.fa", 20, 50, 4),
        //("./data/Synthetic1000.fa", 100, 100, 8),
        //("./data/Synthetic10000.fa", 200, 100, 12),
        //("./data/Synthetic30000.fa", 200, 100, 12),
        //("./data/Synthetic100000.fa", 250, 200, 12),
        //("./data/Pacbio.fa", 200, 100, 12),
    ];
    
    for (filename, bandwidth, cut_off, k) in test_cases {
        let result = bench_mark_one(filename, bandwidth, cut_off, k);
        // time
        println!("Time Taken (lcsk, threaded){:?}", filename, result.0);
        println!("Score (lcsk, threaded){:?}", filename, result.1);
        println!("Mem Used (lcsk, threaded){:?}", filename, result.2);
    }
}

fn bench_mark_one (filename: &str, bandwidth: usize, cut_off: usize, k: usize) -> ((usize, usize), (usize, usize), (usize, usize)) {
    let mut overall_result = ((0, 0), (0, 0), (0, 0));
    let string_vec_vec = get_sequences_from_file(filename);
    for string_vec in string_vec_vec {
        let result = run_test(string_vec, bandwidth, cut_off, k);
        overall_result.0.0 += result.0.0;
        overall_result.0.1 += result.0.1;
        
        overall_result.1.0 += result.1.0;
        overall_result.1.1 += result.1.1;

        overall_result.2.0 += result.2.0;
        overall_result.2.1 += result.2.1;
    }
    overall_result
}

fn run_test (string_vec: Vec<String>, bandwidth: usize, cut_off: usize, k: usize) -> ((usize, usize), (usize, usize), (usize, usize)) {
    let mut lcsk_time_taken = 0;
    let mut lcsk_memory_used = 0;
    let mut lcsk_final_score = 0;

    let mut threaded_time_taken = 0;
    let mut threaded_memory_used = 0;
    let mut threaded_final_score = 0;

    let mut lcsk_aligner = Aligner::new(MATCH_SCORE, MISMATCH_SCORE, -GAP_OPEN_SCORE, &string_vec[0].as_bytes().to_vec());
    let mut threaded_aligner = Aligner::new(MATCH_SCORE, MISMATCH_SCORE, -GAP_OPEN_SCORE, &string_vec[0].as_bytes().to_vec());

    let mut lcsk_all_paths = vec![(0..string_vec[0].len()).collect()];
    let mut lcsk_all_bases = vec![string_vec[0].as_bytes().to_vec()];

    let mut threaded_all_paths = vec![(0..string_vec[0].len()).collect()];
    let mut threaded_all_bases = vec![string_vec[0].as_bytes().to_vec()];

    for (index, string) in string_vec.iter().enumerate().skip(1) {
        let query = &string.as_bytes().to_vec();
        //simd lcskpoa
        let lcsk_graph = lcsk_aligner.graph();
        let now = Instant::now();
        let lcsk_path = lcsk_pipeline (lcsk_graph, query, k, &lcsk_all_paths, &lcsk_all_bases);
        let (obtained_temp_bases, obtained_temp_path, lcsk_score, lcsk_mem_usage) = lcsk_aligner.global_simd_banded(query, &lcsk_path, bandwidth);
        let lcsk_time = now.elapsed().as_micros() as usize;
        lcsk_all_bases.push(obtained_temp_bases);
        lcsk_all_paths.push(obtained_temp_path);
        //threaded simd lcskpoa
        let mut children = vec![];
        let threaded_graph = threaded_aligner.graph();
        let now = Instant::now();
        let (anchors, section_graphs, section_node_trackers, section_queries, section_lcsks) = threaded_lcsk_pipeline (threaded_graph, query, k, &threaded_all_paths, &threaded_all_bases, cut_off);
        let mut threaded_score = 0;
        let mut threaded_mem_usage = 0;
        let mut all_alignment_ops = vec![];
        for anchor_index in 0..anchors.len() + 1 {
            if anchors.len() == 0 {
                threaded_score = lcsk_score;
                // this is just single threaded
                break;
            }
            let mut anchor_start = (0, 0);
            if anchor_index > 0 {
                anchor_start = anchors[anchor_index - 1];
            }
            let section_query = section_queries[anchor_index].clone();
            let section_lcsk = section_lcsks[anchor_index].clone();
            let section_graph = section_graphs[anchor_index].clone();
            let section_node_tracker = section_node_trackers[anchor_index].clone();
            children.push(thread::spawn(move || {
                let mut aligner_child = Aligner::empty(MATCH_SCORE, MISMATCH_SCORE, -GAP_OPEN_SCORE);
                let result = aligner_child.global_simd_banded_threaded_part1 (anchor_start, anchor_index, &section_query, &section_lcsk, bandwidth, section_graph, section_node_tracker);
                result
            }));  
        }
        // concat the alignment
        for _child_index in 0..children.len() {
            let (alignment, mem_usage) = children.remove(0).join().unwrap();
            threaded_score += alignment.score;
            threaded_mem_usage += mem_usage;
            for op in alignment.operations {
                all_alignment_ops.push(op);
            }
        }
        let threaded_time = now.elapsed().as_micros() as usize;
        // for the total score anchors.len() * match score should be added as anchored at matches
        threaded_score = threaded_score + anchors.len() as i32 * MATCH_SCORE;
        println!("threaded score {} {}", threaded_score, lcsk_score);
        let (obtained_temp_bases, obtained_temp_path) = threaded_aligner.global_simd_banded_threaded_part2 (Alignment {score: threaded_score, operations: all_alignment_ops}, query);
        threaded_all_bases.push(obtained_temp_bases);
        threaded_all_paths.push(obtained_temp_path);
        // results we need
        if index == 2 {
            threaded_final_score = threaded_score as usize;
            lcsk_final_score = lcsk_score as usize;

            threaded_time_taken = threaded_time;
            lcsk_time_taken = lcsk_time;

            threaded_memory_used = threaded_mem_usage;
            lcsk_memory_used = lcsk_mem_usage;
        }
    }
    ((lcsk_time_taken, threaded_time_taken), (lcsk_memory_used, threaded_memory_used), (lcsk_final_score, threaded_final_score))
}

fn get_sequences_from_file (filename: &str) -> Vec<Vec<String>> {
    let path = Path::new(filename);
    let file = File::open(&path).unwrap();
    let reader = io::BufReader::new(file);
    
    let mut result: Vec<Vec<String>> = Vec::new();
    let mut temp_vec: Vec<String> = Vec::new();
    
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('>') {
            continue;
        }
        temp_vec.push(line);
        if temp_vec.len() == 3 {
            result.push(temp_vec.clone());
            temp_vec.clear();
        }
    }
    result
}

fn get_random_sequences_from_generator(sequence_length: usize, num_of_sequences: usize, seed: usize) -> Vec<String> {
    let mut rng = StdRng::seed_from_u64(seed as u64);
    //vector to save all the sequences 
    let mut randomvec: Vec<String> = vec![];
    //generate the first sequence of random bases of length sequence_length
    let mut firstseq: Vec<char> = vec![];
    for _ in 0..sequence_length {
        firstseq.push(match rng.gen_range(0..4) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'X'
        });
    }
    //loop for 10 
    for _ in 0..num_of_sequences {
        //clone the sequence
        let mut mutseq = firstseq.clone();
        //mutate the all the bases with 0.05 chance
        for i in 0..mutseq.len() {
            match rng.gen_range(0..20) {
                0 => {
                    mutseq[i] = match rng.gen_range(0..4) {
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        3 => 'T',
                        _ => 'X'
                    }
                },
                _ => {}
            }
        }
        //put indels at location with chance 0.1 
        for i in 0..mutseq.len() {
            let mean_value: f64 = 1.5; //2.0 before
            //get length of the indel geometric distributed mean value 1.5
            let indel_length: usize  = ((1.0 - rng.gen::<f64>()).ln() / (1.00 - (1.00 / mean_value) as f64).ln()).ceil() as usize;
            match rng.gen_range(0..20) {
                //insertion of elements
                0 => {
                    if i + indel_length < mutseq.len() {
                        for _ in 0..indel_length{
                            mutseq.insert(i + 1, mutseq[i]);
                        }
                    }
                },
                //deletion of elements
                1 => {
                    if i + indel_length < mutseq.len() {
                        for _ in 0..indel_length{
                            mutseq.remove(i);
                        }
                    }
                }
                _ => {}
            }
        }
        //insert to vector
        randomvec.push(mutseq.iter().collect::<String>());
    }
    randomvec
}