#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
use poa::*;
use std::simd::{i32x4, Simd};

fn main() {
    // TEST SIMD Add two arrays of 4 elements each in parallel using SIMD
    // Create two arrays
    let a = [1, 2, 3, 4];
    let b = [5, 6, 7, 8];

    // Load the arrays into SIMD vectors (f32x4 can hold 4 floats in parallel)
    let simd_a = i32x4::from_array(a);
    let simd_b = i32x4::from_array(b);

    // Perform element-wise addition using SIMD
    let simd_result = simd_a + simd_b;

    // Extract the result back into a regular array
    let result: [i32; 4] = simd_result.to_array();

    // Print the result
    println!("Result: {:?}", result);
    
    // TEST POA
    let seqs = ["ACT", "AGT", "ACC", "ACT"];
    let mut seqs_bytes = vec![];
    for seq in seqs.iter() {
        seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
    }
    let mut aligner = Aligner::new(2, -2, -2, &seqs_bytes[0].to_vec(), 0, 0, 10 as i32);
    for seq in seqs_bytes.iter().skip(1) {
        aligner.global(seq).add_to_graph();
    }
    let consensus = aligner.consensus();
    println!("{:?}", consensus)
}