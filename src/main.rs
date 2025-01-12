#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
mod pairwise;
use pairwise::*;
use poa::*;
use std::simd::{u16x4, Simd};

fn main() {
    let seq_x = vec![65, 65, 84, 65, 65, 84, 65, 65];
    let seq_y = vec![65, 65, 84, 65, 65, 84, 65, 65];
    let match_score = 1;
    let mismatch_score = -1;
    let gap_open_score = 2;
    let gap_extend_score = 1;
    pairwise_simd(&seq_x, &seq_y, match_score, mismatch_score, gap_open_score, gap_extend_score);
    //fake_pairwise_simd(&seq_x, &seq_y, match_score, mismatch_score, gap_open_score, gap_extend_score);
    //pairwise(&seq_x, &seq_y, match_score, mismatch_score, -gap_open_score, -gap_extend_score, 10);
}