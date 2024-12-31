
use std::cmp;
const MIN_SCORE: isize = -858_993_459; // negative infinity
use std::simd::{i32x4, Simd};

#[derive(Clone)]
struct PairwiseMatrixCell {
    match_score: isize,
    del_score: isize,
    ins_score: isize,
    back: char,
}

pub fn pairwise_simd (seq_x: &Vec<u8>, seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32) -> (Vec<u8>, isize) {
    // variables to save results
    let mut align_vec: Vec<u8> = Vec::new();

    // make one matrix 
    let mut pair_wise_matrix: Vec<Vec<PairwiseMatrixCell>> = vec![vec![PairwiseMatrixCell { match_score: (0), del_score: (0), ins_score: (0), back: ('X') }; seq_y.len() + 1]; seq_x.len() + 1];

    // fill out the first row and colomn
    let temp_value = gap_open_score as isize + gap_extend_score as isize;
    pair_wise_matrix[0][1] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('d') };
    pair_wise_matrix[1][0] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('i') };
    for i in 2..seq_y.len() + 1 {
        let temp_value = pair_wise_matrix[0][i - 1].del_score + gap_extend_score as isize;
        pair_wise_matrix[0][i] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('d') };
    }
    for i in 2..seq_x.len() + 1 {
        let temp_value = pair_wise_matrix[i - 1][0].ins_score + gap_extend_score as isize;
        pair_wise_matrix[i][0] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('i') };
    }
    // calculations
    // filling out score matrices and back matrix
    for i in 1..seq_x.len() + 1 {
        for j in 1..seq_y.len() + 1 {
            // fill del matrix 
            // get j - 1 score from same matrix with gap extend
            let temp_del_score = pair_wise_matrix[i][j - 1].del_score + gap_extend_score as isize;
            // get j - 1 score from match matrix with gap open penalty
            let temp_match_score = pair_wise_matrix[i][j - 1].match_score + gap_open_score as isize + gap_extend_score as isize;
            // insert the max
            pair_wise_matrix[i][j].del_score = cmp::max(temp_del_score, temp_match_score);

            // fill ins matrix
            // get i - 1 score from the same matrix
            let temp_ins_score = pair_wise_matrix[i - 1][j].ins_score + gap_extend_score as isize;
            // get i - 1 score from the match matrix with gap open penalty
            let temp_match_score = pair_wise_matrix[i - 1][j].match_score + gap_open_score as isize + gap_extend_score as isize;
            // insert the max
            pair_wise_matrix[i][j].ins_score = cmp::max(temp_ins_score, temp_match_score);

            // fill match matrix
            // get the i,j from the insertion matrix
            let temp_ins_score = pair_wise_matrix[i][j].ins_score;
            // get the i,j from the deletion matrix
            let temp_del_score = pair_wise_matrix[i][j].del_score;
            // get the match from i-1,j-1 from match matrix with match score or mismatch score
            let temp_match_score;
            if seq_x[i - 1] == seq_y[j - 1] {
                temp_match_score = pair_wise_matrix[i - 1][j - 1].match_score + match_score as isize;
                pair_wise_matrix[i][j].back = 'm';
            }
            else {
                temp_match_score = pair_wise_matrix[i - 1][j - 1].match_score + mismatch_score as isize;
                pair_wise_matrix[i][j].back = 's';
            }
            // insert the max
            pair_wise_matrix[i][j].match_score = cmp::max(temp_match_score, cmp::max(temp_ins_score, temp_del_score));
            if (temp_match_score >= temp_ins_score) && (temp_match_score >= temp_del_score) {
                // already allocated
            }
            else if temp_ins_score > temp_del_score {
                pair_wise_matrix[i][j].back = 'i';
            }
            else {
                pair_wise_matrix[i][j].back = 'd';
            }
        }
    }
    // back tracing using back matrix and filling out align_vec
    let mut i = seq_x.len();
    let mut j = seq_y.len();
    let score = pair_wise_matrix[i][j].match_score;
    let mut break_on_next = false;
    // print the score matrix
    println!("Score matrix");
    for i_index in 0..seq_x.len() {
        for j_index in 0..seq_y.len() {
            print!(" {}", pair_wise_matrix[i_index][j_index].match_score);
        }
        println!("");
    }
    // print the del matrix
    println!("Del matrix");
    for i_index in 0..seq_x.len() {
        for j_index in 0..seq_y.len() {
            print!(" {}", pair_wise_matrix[i_index][j_index].del_score);
        }
        println!("");
    }
    // print the ins matrix
    println!("Ins matrix");
    for i_index in 0..seq_x.len() {
        for j_index in 0..seq_y.len() {
            print!(" {}", pair_wise_matrix[i_index][j_index].ins_score);
        }
        println!("");
    }
    loop {
        match pair_wise_matrix[i][j].back {
            'i' => {
                i = i - 1;
                align_vec.push('i' as u8);
            },
            'm' => {
                i = i - 1;
                j = j - 1;
                align_vec.push('m' as u8);
            },
            's' => {
                i = i - 1;
                j = j - 1;
                align_vec.push('s' as u8);
            }
            'd' => {
                j = j - 1;
                align_vec.push('d' as u8);
            },
            _ => (),
        }
        if break_on_next {
            break;
        }
        if i == 0 && j == 0 {
            break_on_next = true;
        }
    }
    (align_vec.into_iter().rev().collect(), score)
}

pub fn pairwise (seq_x: &Vec<u8>, seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32, band_size: usize) -> (Vec<u8>, isize) {
    // variables to save results
    let mut align_vec: Vec<u8> = Vec::new();

    // make one matrix 
    let mut pair_wise_matrix: Vec<Vec<PairwiseMatrixCell>> = vec![vec![PairwiseMatrixCell { match_score: (0), del_score: (0), ins_score: (0), back: ('X') }; seq_y.len() + 1]; seq_x.len() + 1];

    // fill out the first row and colomn
    let temp_value = gap_open_score as isize + gap_extend_score as isize;
    pair_wise_matrix[0][1] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('d') };
    pair_wise_matrix[1][0] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('i') };
    for i in 2..seq_y.len() + 1 {
        let temp_value = pair_wise_matrix[0][i - 1].del_score + gap_extend_score as isize;
        pair_wise_matrix[0][i] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('d') };
    }
    for i in 2..seq_x.len() + 1 {
        let temp_value = pair_wise_matrix[i - 1][0].ins_score + gap_extend_score as isize;
        pair_wise_matrix[i][0] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('i') };
    }
    // calculations
    // filling out score matrices and back matrix
    let mut max_scored_position = 0;
    let mut max_score;
    let mut start = 0;
    let mut end= seq_x.len();
    for i in 1..seq_x.len() + 1 {
        max_score = MIN_SCORE;
        if band_size > 0 {
            if max_scored_position < band_size {
                start = 0;
            }
            else {
                start = max_scored_position - band_size;
            }
            end = max_scored_position + band_size;
            // start at 0 initially
            if i < 20 {
                start = 0;
            }
            // end at end at end :D
            if i > seq_x.len() - 20 {
                end = seq_x.len();
            }
        }
        
        for j in 1..seq_y.len() + 1 {
            if j < start && (band_size > 0) && i > 1 {
                continue;
            }
            if j > end && (band_size > 0) && i > 1{
                break;
            }
            // fill del matrix 
            // get j - 1 score from same matrix with gap extend
            let temp_del_score = pair_wise_matrix[i][j - 1].del_score + gap_extend_score as isize;
            // get j - 1 score from match matrix with gap open penalty
            let temp_match_score = pair_wise_matrix[i][j - 1].match_score + gap_open_score as isize + gap_extend_score as isize;
            // insert the max
            pair_wise_matrix[i][j].del_score = cmp::max(temp_del_score, temp_match_score);

            // fill ins matrix
            // get i - 1 score from the same matrix
            let temp_ins_score = pair_wise_matrix[i - 1][j].ins_score + gap_extend_score as isize;
            // get i - 1 score from the match matrix with gap open penalty
            let temp_match_score = pair_wise_matrix[i - 1][j].match_score + gap_open_score as isize + gap_extend_score as isize;
            // insert the max
            pair_wise_matrix[i][j].ins_score = cmp::max(temp_ins_score, temp_match_score);

            // fill match matrix
            // get the i,j from the insertion matrix
            let temp_ins_score = pair_wise_matrix[i][j].ins_score;
            // get the i,j from the deletion matrix
            let temp_del_score = pair_wise_matrix[i][j].del_score;
            // get the match from i-1,j-1 from match matrix with match score or mismatch score
            let temp_match_score;
            if seq_x[i - 1] == seq_y[j - 1] {
                temp_match_score = pair_wise_matrix[i - 1][j - 1].match_score + match_score as isize;
                pair_wise_matrix[i][j].back = 'm';
            }
            else {
                temp_match_score = pair_wise_matrix[i - 1][j - 1].match_score + mismatch_score as isize;
                pair_wise_matrix[i][j].back = 's';
            }
            // insert the max
            pair_wise_matrix[i][j].match_score = cmp::max(temp_match_score, cmp::max(temp_ins_score, temp_del_score));
            if max_score < pair_wise_matrix[i][j].match_score {
                max_score = pair_wise_matrix[i][j].match_score;
                max_scored_position = j;
            }
            if (temp_match_score >= temp_ins_score) && (temp_match_score >= temp_del_score) {
                // already allocated
            }
            else if temp_ins_score > temp_del_score {
                pair_wise_matrix[i][j].back = 'i';
            }
            else {
                pair_wise_matrix[i][j].back = 'd';
            }
        }
    }
    // back tracing using back matrix and filling out align_vec
    let mut i = seq_x.len();
    let mut j = seq_y.len();
    let score = pair_wise_matrix[i][j].match_score;
    let mut break_on_next = false;
    // print the score matrix
    println!("Score matrix");
    for i_index in 0..seq_x.len() {
        for j_index in 0..seq_y.len() {
            print!(" {}", pair_wise_matrix[i_index][j_index].match_score);
        }
        println!("");
    }
    // print the del matrix
    println!("Del matrix");
    for i_index in 0..seq_x.len() {
        for j_index in 0..seq_y.len() {
            print!(" {}", pair_wise_matrix[i_index][j_index].del_score);
        }
        println!("");
    }
    // print the ins matrix
    println!("Ins matrix");
    for i_index in 0..seq_x.len() {
        for j_index in 0..seq_y.len() {
            print!(" {}", pair_wise_matrix[i_index][j_index].ins_score);
        }
        println!("");
    }
    loop {
        match pair_wise_matrix[i][j].back {
            'i' => {
                i = i - 1;
                align_vec.push('i' as u8);
            },
            'm' => {
                i = i - 1;
                j = j - 1;
                align_vec.push('m' as u8);
            },
            's' => {
                i = i - 1;
                j = j - 1;
                align_vec.push('s' as u8);
            }
            'd' => {
                j = j - 1;
                align_vec.push('d' as u8);
            },
            _ => (),
        }
        if break_on_next {
            break;
        }
        if i == 0 && j == 0 {
            break_on_next = true;
        }
    }
    (align_vec.into_iter().rev().collect(), score)
}