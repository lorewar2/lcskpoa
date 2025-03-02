#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
mod poa_old;
mod lcsk;
use poa_old::Aligner as Aligner2;
use petgraph::{dot::Dot, Direction::Outgoing};
use poa::*;
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::thread;
use petgraph::Incoming;
use petgraph::graph::NodeIndex;
mod bit_tree;
use crate::lcsk::{lcsk_pipeline, threaded_lcsk_pipeline};

fn main() {
    for seed in 0..100 {
        println!("seed {}", seed);
        let seqs = get_random_sequences_from_generator(100, 8, seed);
        let match_score = 2;
        let mismatch_score = -2;
        let gap_open_score = 2;
        let band_size = 20;
        let kmer_size = 8;
        let cut_limit = 50;
        let mut children = vec![];
        // full aligner
        let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec());
        let mut old_aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec());
        //let mut old_aligner = Aligner2::new(match_score, mismatch_score, -gap_open_score, &seqs[0].as_bytes().to_vec(), i32::MIN, i32::MIN, band_size);
        // you have to keep the original graph otherwise these things change
        let mut all_paths = vec![(0..seqs[0].len()).collect()];
        let mut all_bases = vec![seqs[0].as_bytes().to_vec()];
        for index in 1..seqs.len() - 1 {
            println!("Running {} seq {:?}", index, seqs[index]);
            let query = &seqs[index].as_bytes().to_vec();
            let output_graph = aligner.graph();
            //println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
            let lcsk_path = lcsk_pipeline (output_graph, query, kmer_size, &all_paths, &all_bases);
            //println!("lcsk path {:?}", lcsk_path);
            let (anchors, section_graphs, section_node_trackers, section_queries, section_lcsks) = threaded_lcsk_pipeline (output_graph, query, kmer_size, &all_paths, &all_bases, cut_limit);
            let mut total_score = 0;
            let mut all_alignment_ops = vec![];
            for anchor_index in 0..anchors.len() + 1 {
                if anchors.len() == 0 {
                    break;
                }
                let section_query = section_queries[anchor_index].clone();
                let section_lcsk = section_lcsks[anchor_index].clone();
                let section_graph = section_graphs[anchor_index].clone();
                let section_node_tracker = section_node_trackers[anchor_index].clone();
                //println!(" {} {:?} {} {}", anchors.len(), section_query, section_lcsks.len(), section_graphs.len());
                for base in &section_query {
                    //print!("{}", *base as char);
                }
                //println!("\nnumber of nodes per section {}", section_graph.node_count());
                //println!("query len {}, graph len {}", section_query.len(), section_graph.node_count());
                //println!("{:?}", Dot::new(&section_graph.map(|_, n| (*n) as char, |_, e| *e)));
                children.push(thread::spawn(move || {
                    let mut aligner2 = Aligner::empty(2, -2, -2, band_size);
                    let alignment = aligner2.global_simd_banded_threaded_part1(&section_query, &section_lcsk, band_size as usize, section_graph, section_node_tracker);
                    alignment
                }));  
            }
            // concat the alignment
            for child_index in 0..children.len() {
                let alignment = children.remove(0).join().unwrap();
                println!("current_score {}", alignment.score);
                //println!("alignment {:?}", alignment.operations);
                total_score += alignment.score;
                for op in alignment.operations {
                    if op == AlignmentOperation::Match(None) && child_index > 0 {
                        // this should only have one parent as this is an anchor
                        let incoming_nodes: Vec<NodeIndex<usize>> = output_graph.neighbors_directed(NodeIndex::new(anchors[child_index - 1]), Outgoing).collect();
                        assert!(incoming_nodes.len() == 1);
                        all_alignment_ops.push(AlignmentOperation::Match(Some((anchors[child_index - 1], incoming_nodes[0].index()))));
                    }
                    else {
                        all_alignment_ops.push(op);
                    }
                }
            }
            // using the alignment add to the graph and get stuff
            //println!("section combined assignment {:?}", all_alignment_ops);
            // for the total score anchors.len() * match score should be added as anchored at matches
            println!("threaded score {}", total_score + anchors.len() as i32 * match_score);
            let (obtained_temp_bases, obtained_temp_path) = aligner.global_simd_banded_threaded_part2(Alignment {score: total_score, operations: all_alignment_ops}, query);
            //let (obtained_temp_bases, obtained_temp_path) = old_aligner.global_simd_banded(query, &lcsk_path, band_size as usize);
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