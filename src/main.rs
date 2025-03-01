#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
mod poa_old;
mod lcsk;
use poa_old::Aligner as Aligner2;
use petgraph::dot::Dot;
use poa::*;
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::thread;
mod bit_tree;
use crate::lcsk::{lcsk_pipeline, threaded_lcsk_pipeline};

fn main() {
    for seed in 0..1 {
        println!("seed {}", seed);
        let seqs = get_random_sequences_from_generator(1000, 8, seed);
        let match_score = 2;
        let mismatch_score = -2;
        let gap_open_score = 2;
        let band_size = 100;
        let kmer_size = 8;
        let cut_limit = 200;
        let mut children = vec![];
        let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec());
        //let mut old_aligner = Aligner2::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec(), i32::MIN, i32::MIN, band_size);
        // you have to keep the original graph otherwise these things change
        let mut all_paths = vec![(0..seqs[0].len()).collect()];
        let mut all_bases = vec![seqs[0].as_bytes().to_vec()];
        for index in 1..seqs.len() - 1 {
            println!("Running {} seq {:?}", index, seqs[index]);
            let query = &seqs[index].as_bytes().to_vec();
            let output_graph = aligner.graph();
            println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
            let lcsk_path = lcsk_pipeline (output_graph, query, kmer_size, &all_paths, &all_bases);
            println!("lcsk path {:?}", lcsk_path);
            let (anchors, section_graphs, _node_tracker, section_queries, section_lcsks) = threaded_lcsk_pipeline (output_graph, query, kmer_size, &all_paths, &all_bases, cut_limit);
            
            let mut total_score = 0;
            for anchor_index in 0..anchors.len() + 1 {
                if anchors.len() == 0 {
                    break;
                }
                let section_query = section_queries[anchor_index].clone();
                let section_lcsk = section_lcsks[anchor_index].clone();
                let section_graph = section_graphs[anchor_index].clone();
                //println!(" {} {:?} {} {}", anchors.len(), section_query, section_lcsks.len(), section_graphs.len());
                for base in &section_query {
                    print!("{}", *base as char);
                }
                println!("\nnumber of nodes per section {}", section_graph.node_count());
                //println!("query len {}, graph len {}", section_query.len(), section_graph.node_count());
                println!("{:?}", Dot::new(&section_graph.map(|_, n| (*n) as char, |_, e| *e)));
                
                children.push(thread::spawn(move || {
                    let mut aligner = Aligner::empty(2, -2, -2, band_size);
                    let score = aligner.global_simd_banded_threaded(&section_query, &section_lcsk, band_size as usize, section_graph);
                    score
                }));
                for _child_index in 0..children.len() {
                    let result = children.pop().unwrap().join().unwrap();
                    println!("current_score {}", result);
                    total_score += result;
                }
            }
            println!("threaded score {}", total_score);
            let (obtained_temp_bases, obtained_temp_path) = aligner.global_simd_banded(query, &lcsk_path, band_size as usize);
            //let (obtained_temp_bases, obtained_temp_path) = aligner.global_simd(query);
            //old_aligner.custom(query).add_to_graph();
            all_bases.push(obtained_temp_bases);
            all_paths.push(obtained_temp_path);
            
        }
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