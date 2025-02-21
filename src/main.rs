#![feature(portable_simd)]
//#![allow(dead_code)]
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
use crate::lcsk::{find_sequence_in_graph, find_kmer_matches, lcskpp_graph};

fn main() {
    let seqs = get_random_sequences_from_generator(12, 8, 0);
    let match_score = 2;
    let mismatch_score = -2;
    let gap_open_score = 2;
    // get the graph first
    let band_size = 250;
    let seed = 1;
    let kmer_size = 3;
    let mut all_paths = vec![];
    let mut all_sequences = vec![];
    //let seqs = get_random_sequences_from_generator(1000, 3, seed);
    println!("seqs1 {}", seqs[4]);
    //let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec());
    let mut old_aligner = Aligner2::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec(), i32::MIN, i32::MIN, &band_size);
    for index in 1..seqs.len() - 1 {
        let output_graph = old_aligner.graph();
        println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
        println!("{}", seqs[index]);
        let mut topo = Topo::new(&output_graph);
        let mut topo_indices = vec![];
        let mut topo_map = HashMap::new();
        let mut incrementing_index: usize = 0;
        while let Some(node) = topo.next(&output_graph) {
            print!("{} {}, ", node.index(), output_graph.raw_nodes()[node.index()].weight as char);
            topo_indices.push(node.index());
            
            topo_map.insert(node.index(), incrementing_index);
            incrementing_index += 1;
        }
        println!("");
        let mut error_index = 0;
        loop {
            let (error_occured, temp_path, temp_sequence) = find_sequence_in_graph(seqs[index - 1].as_bytes().to_vec().clone(), &output_graph, &topo_indices, &topo_map, error_index);
            if error_index > 10 {
                break;
            }
            if !error_occured {
                println!("temp path {:?}", temp_path);
                println!("temp seq {:?}", temp_sequence);
                all_paths.push(temp_path);
                all_sequences.push(temp_sequence);
                println!("NO error occured");
                break;
            }
            error_index += 1;   
        }
        println!("all paths vec len {}", all_paths.len());
        let query = &seqs[index].as_bytes().to_vec();
        println!("query {:?}", query);
        println!("path of last {:?}", all_paths.last());
        let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = find_kmer_matches(&query, &all_sequences, &all_paths, kmer_size);
        println!("kmers {:?} len {}", kmer_pos_vec, kmer_pos_vec.len());
        for index_2 in 0..kmer_pos_vec.len() {
            println!("node index {} original kmer {:?} {:?} {:?} {:?}", topo_indices[kmer_pos_vec[index_2].1 as usize], kmer_pos_vec[index_2], kmer_path_vec[index_2], kmers_previous_node_in_paths[index_2], kmer_graph_path[index_2]);
        }
        //let (lcsk_path, _lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), kmer_size, kmer_graph_path, &topo_indices);
        //println!("{:?}", lcsk_path);
        let now = Instant::now();
        //aligner.global_simd(query);

        old_aligner.global(query).add_to_graph();

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