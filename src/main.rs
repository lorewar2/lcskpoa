#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
use poa::*;
use std::time::Instant;
use petgraph::dot::Dot;
use rand::{Rng, SeedableRng, rngs::StdRng};

fn main() {
    let match_score = 1;
    let mismatch_score = -1;
    let gap_open_score = 2;
    for seed in 0..1 {
        let seqs = get_random_sequences_from_generator(10000, 2, seed);
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
        }
        let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs_bytes[0]);
        let now = Instant::now();
        for seq in seqs_bytes.iter().skip(1) {
            aligner.global(seq).add_to_graph();
        }
        let graph = aligner.graph();
        //println!("{:?}", Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e)));
        let time = now.elapsed().as_micros() as usize;
        println!("Completed normal poa elapsed time {}μs", time);
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
        }
        let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs_bytes[0]);
        let now = Instant::now();
        for seq in seqs_bytes.iter().skip(1) {
            aligner.global_simd(seq);
        }
        let graph = aligner.graph();
        //println!("{:?}", Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e)));
        let time = now.elapsed().as_micros() as usize;
        println!("Completed simd poa elapsed time {}μs", time);
    }
    
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
    //randomvec.push(firstseq.iter().collect::<String>());
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
        //println!("{:?}", mutseq.iter().collect::<String>());
        //insert to vector
        randomvec.push(mutseq.iter().collect::<String>());
    }
    randomvec
}