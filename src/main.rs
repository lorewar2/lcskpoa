#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
mod poa_old;
use poa_old::Aligner as Aligner2;
use poa::*;
use std::time::Instant;
use rand::{Rng, SeedableRng, rngs::StdRng};
mod lcsk;
mod bit_tree;
use std::collections::HashMap;
use petgraph::visit::Topo;
use crate::lcsk::{find_sequence_in_graph, better_find_kmer_matches, lcskpp_graph, anchoring_lcsk_path_for_threading};

fn main() {
    let seqs = vec![
        "AAAAT",
        "AAAAA",
        "AAAACGT"
    ];
    let match_score = 2;
    let mismatch_score = -2;
    let gap_open_score = 2;
    // get the graph first
    let band_size = 250;
    let seed = 1;
    let kmer_size = 12;
    let mut all_paths = vec![];
    let mut all_sequences = vec![];
    //let seqs = get_random_sequences_from_generator(1000, 3, seed);
    //println!("seqs1 {}", seqs[0]);
    let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec());
    let mut old_aligner = Aligner2::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec(), i32::MIN, i32::MIN, &band_size);
    for index in 1..seqs.len() {
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
        //println!("{:?}", lcsk_path);
        let now = Instant::now();
        //aligner.global_simd(query);
        if index != 20 {
            //old_aligner.global_banded(query, &lcsk_path, band_size as usize).add_to_graph();
            aligner.global_simd(query);
            old_aligner.global(query).add_to_graph();
            println!("{}", old_aligner.alignment().score);
            //aligner.global_simd_banded(query, &lcsk_path, band_size as usize);
        }
        else {
            //old_aligner.global(query).add_to_graph();
            //aligner.global(query).add_to_graph();
        }
        let time = now.elapsed().as_micros() as usize;
        println!("Completed kmer time elapsed time {}μs", time);
    }
}

fn read_from_fa_give_vec (path: String) -> Vec<Vec<String>> {
    let mut string_vec = vec![];

    string_vec
}