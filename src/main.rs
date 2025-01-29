#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
use poa::*;
use std::time::Instant;
use rand::{Rng, SeedableRng, rngs::StdRng};
mod lcsk;
mod bit_tree;
use std::collections::HashMap;
use petgraph::visit::Topo;
use crate::lcsk::{find_sequence_in_graph, better_find_kmer_matches, lcskpp_graph, anchoring_lcsk_path_for_threading};

fn main() {
    let match_score = 1;
    let mismatch_score = -1;
    let gap_open_score = 2;
    // get the graph first
    
    let seed = 0;
    let kmer_size = 4;
    let mut all_paths = vec![];
    let mut all_sequences = vec![];
    let seqs = get_random_sequences_from_generator(400_000, 10, seed);
    let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec());
    for index in 1..seqs.len() {
        let now = Instant::now();
        let output_graph = aligner.graph();
        let mut topo = Topo::new(&output_graph);
        let mut topo_indices = vec![];
        let mut topo_map = HashMap::new();
        let mut incrementing_index: usize = 0;
        while let Some(node) = topo.next(&output_graph) {
            topo_indices.push(node.index());
            topo_map.insert(node.index(), incrementing_index);
            incrementing_index += 1;
        }
        let mut error_index = 0;
        loop {
            let (error_occured, temp_path, temp_sequence) = find_sequence_in_graph (seqs[index - 1].as_bytes().to_vec().clone(), &output_graph, &topo_indices, &topo_map, error_index);
            if error_index > 10 {
                break;
            }
            if !error_occured {
                all_paths.push(temp_path);
                all_sequences.push(temp_sequence);
                break;
            }
            error_index += 1;   
        }
        let query = &seqs[index].as_bytes().to_vec();
        let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&query, &all_sequences, &all_paths, kmer_size);
        let (lcsk_path, _lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), kmer_size, kmer_graph_path, &topo_indices);
        let time = now.elapsed().as_micros() as usize;
        println!("Completed kmer time elapsed time {}μs", time);
        //aligner.global_simd(query);
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