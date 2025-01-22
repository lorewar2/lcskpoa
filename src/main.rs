#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
use pairwise::*;
use poa::*;
use std::simd::{u16x4, Simd};
use std::time::Instant;

fn main() {
    let seq_x = vec![65, 65, 84, 65, 65, 84, 65, 65];
    let seq_y = vec![65, 65, 84, 65, 65, 84, 65, 65];
    let match_score = 1;
    let mismatch_score = -1;
    let gap_open_score = 2;
    let gap_extend_score = 1;
    //pairwise_simd_without_extend (&seq_x, &seq_y, match_score, mismatch_score, gap_open_score, gap_extend_score);
    //pairwise_without_extend(&seq_x, &seq_y, match_score, mismatch_score, -gap_open_score, -gap_extend_score);
    //fake_pairwise_simd(&seq_x, &seq_y, match_score, mismatch_score, gap_open_score, gap_extend_score);
    //pairwise(&seq_x, &seq_y, match_score, mismatch_score, -gap_open_score, -gap_extend_score, 10);
    // nw for now
    // test poa simple here
    let seqs = vec!["TCTTTCTC".to_string(), "TTCTTTCCC".to_string()];
    let mut seqs_bytes = vec![];
    for seq in seqs.iter() {
        seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
    }
    let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs_bytes[0]);
    let now = Instant::now();
    for seq in seqs_bytes.iter().skip(1) {
        aligner.global(seq);
    }
    let time = now.elapsed().as_micros() as usize;
    println!("Completed normal poa elapsed time {}μs", time);

    let seqs = vec!["TCTTTCTC".to_string(), "TTCTTTCCC".to_string()];
    let mut seqs_bytes = vec![];
    for seq in seqs.iter() {
        seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
    }
    let mut aligner = Aligner::new(match_score, mismatch_score, -gap_open_score, &seqs_bytes[0]);
    let now = Instant::now();
    for seq in seqs_bytes.iter().skip(1) {
        aligner.global_simd(seq);
    }
    
    let time = now.elapsed().as_micros() as usize;
    println!("Completed simd poa elapsed time {}μs", time);
}