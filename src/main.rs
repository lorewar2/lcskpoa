#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
mod poa_old;
use poa_old::Aligner as Aligner2;
use poa::*;
use petgraph::dot::Dot;
use std::time::Instant;
use rand::{Rng, SeedableRng, rngs::StdRng};
mod lcsk;
mod bit_tree;
use std::collections::HashMap;
use petgraph::visit::Topo;
use crate::lcsk::{find_kmer_matches, lcskpp_graph};

fn main() {
    let seqs = get_random_sequences_from_generator(12, 8, 0);
    let match_score = 2;
    let mismatch_score = -2;
    let gap_open_score = 2;
    let band_size = 250;
    let seed = 1;
    let kmer_size = 3;
    
    //let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec());
    let mut old_aligner = Aligner2::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec(), i32::MIN, i32::MIN, band_size);
    // you have to keep the original graph otherwise these things change
    let mut temp_vec: Vec<usize> = (0..seqs[0].len()).collect();
    let mut all_paths = vec![temp_vec];
    let mut all_bases = vec![seqs[0].as_bytes().to_vec()];
    for index in 1..seqs.len() - 1 {
        let output_graph = old_aligner.graph();
        println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
        println!("{:?}", seqs[index]);
        let mut topo = Topo::new(&output_graph);
        let mut index_topo_index_value_graph_index_vec = vec![];
        let mut key_graph_index_value_topo_index_map = HashMap::new();
        let mut incrementing_index: usize = 0;
        let mut temp_paths_converted = vec![];
        // convert the all paths to topo indices for the current iteration of graph
        while let Some(node) = topo.next(&output_graph) {
            index_topo_index_value_graph_index_vec.push(node.index());
            key_graph_index_value_topo_index_map.insert(node.index(), incrementing_index);
            incrementing_index += 1;
        }
        for path in &all_paths {
            let mut temp_path = vec![];
            for node_id in path {
                temp_path.push(key_graph_index_value_topo_index_map[node_id]);
            }
            temp_paths_converted.push(temp_path);
        }
        let query = &seqs[index].as_bytes().to_vec();
        for index_2 in 0..all_bases.len() {
            let mut converted = vec![];
            for entry in &temp_paths_converted[index_2] {
                converted.push(index_topo_index_value_graph_index_vec[*entry as usize]);
            }
            println!("{} {:?} {:?}", index_2, all_bases[index_2], converted);
        }
        let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = find_kmer_matches(&query, &all_bases, &temp_paths_converted, kmer_size);
        for index_1 in 0..kmer_path_vec.len() {
            let mut converted = vec![];
            for entry in &kmer_graph_path[index_1] {
                converted.push(index_topo_index_value_graph_index_vec[*entry as usize]);
            }
            println!("{:?} {:?} {:?}",kmer_pos_vec[index_1], kmer_path_vec[index_1], converted);
        }
        let (lcsk_path, _lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, temp_paths_converted.len(), kmer_size, kmer_graph_path, &index_topo_index_value_graph_index_vec);
        println!("{:?}", lcsk_path);
        let now = Instant::now();
        //aligner.global_simd(query);
        // we get the added path sequence and graph node indices when adding to graph
        let (obtained_temp_bases, obtained_temp_path) = old_aligner.custom(query).add_to_graph();
        println!("Obtained bases at end {:?}", obtained_temp_bases);
        println!("Obtained path at end {:?}", obtained_temp_path);
        all_bases.push(obtained_temp_bases);
        all_paths.push(obtained_temp_path);
        let time = now.elapsed().as_micros() as usize;
        println!("Completed kmer time elapsed time {}μs", time);
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