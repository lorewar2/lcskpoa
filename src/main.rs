#![feature(portable_simd)]
#![allow(dead_code)]
mod poa;
use poa::*;
use std::simd::{f32x4, Simd};

fn main() {
    // TEST SIMD Add two arrays of 4 elements each in parallel using SIMD
    // Create two arrays
    let a = [1.0, 2.0, 3.0, 4.0];
    let b = [5.0, 6.0, 7.0, 8.0];

    // Load the arrays into SIMD vectors (f32x4 can hold 4 floats in parallel)
    let simd_a = f32x4::from_array(a);
    let simd_b = f32x4::from_array(b);

    // Perform element-wise addition using SIMD
    let simd_result = simd_a + simd_b;

    // Extract the result back into a regular array
    let result: [f32; 4] = simd_result.to_array();

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