#![feature(portable_simd)]

use std::simd::{f32x4, Simd};

fn main() {
    // Example: Add two arrays of 4 elements each in parallel using SIMD

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
}