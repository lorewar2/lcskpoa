
use std::cmp;
const MIN_SCORE: isize = -858_993_459; // negative infinity
use std::simd::{i16x8, Simd};
use std::simd::cmp::SimdOrd;
use std::simd::cmp::SimdPartialOrd;
use std::collections::HashMap;

#[derive(Clone)]
struct PairwiseMatrixCell {
    match_score: isize,
    del_score: isize,
    ins_score: isize,
    back: char,
}

pub fn profile_query (seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32) -> Vec<Vec<i16x8>> {
    let num_seq_vec = seq_y.len() / 8;
    let mut MM_simd = vec![];
    // make 4 vectors for query
    let mut A_simd: Vec<i16x8> = vec![];
    let mut C_simd: Vec<i16x8> = vec![];
    let mut G_simd: Vec<i16x8> = vec![];
    let mut T_simd: Vec<i16x8> = vec![];
    // go through the query and populate the entries
    for index_simd in 0..num_seq_vec {
        let mut temp_A = vec![];
        let mut temp_C = vec![];
        let mut temp_G = vec![];
        let mut temp_T = vec![];
        let simd_seq = &seq_y[(index_simd * 8)..((index_simd + 1) * 8)];
        for base in simd_seq {
            if *base == 65 {
                temp_A.push(match_score as i16);
            }
            else {
                temp_A.push(mismatch_score as i16);
            }
            if *base == 67 {
                temp_C.push(match_score as i16);
            }
            else {
                temp_C.push(mismatch_score as i16);
            }
            if *base == 71 {
                temp_G.push(match_score as i16);
            }
            else {
                temp_G.push(mismatch_score as i16);
            }
            if *base == 84 {
                temp_T.push(match_score as i16);
            }
            else {
                temp_T.push(mismatch_score as i16);
            }
        }
        A_simd.push(i16x8::from_array(temp_A[0..8].try_into().expect("")));
        C_simd.push(i16x8::from_array(temp_C[0..8].try_into().expect("")));
        G_simd.push(i16x8::from_array(temp_G[0..8].try_into().expect("")));
        T_simd.push(i16x8::from_array(temp_T[0..8].try_into().expect("")));
        //println!("{:?}", simd_seq);
    }
    MM_simd.push(A_simd);
    MM_simd.push(C_simd);
    MM_simd.push(G_simd);
    MM_simd.push(T_simd);
    MM_simd
}

pub fn pairwise_simd (seq_x: &Vec<u8>, seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32) {
    println!("real");
    // initialize hash map for profiling 
    let mut hash_table = HashMap::new();
    hash_table.insert(65, 0);
    hash_table.insert(67, 1);
    hash_table.insert(71, 2);
    hash_table.insert(84, 3);
    // make the MM table 
    let MM_simd_full = profile_query(seq_y, match_score, mismatch_score);

    // filling out score matrices and back matrix
    // the number of vectors for query SIMD 8 for now
    let num_seq_vec = seq_y.len() / 8;
    // make vectors of 8 of size num_seq_vec
    let gap_open_score = gap_open_score as i16;
    let gap_extend_score = gap_extend_score as i16;
    let gap_open_8 = i16x8::from_array([gap_open_score, gap_open_score, gap_open_score, gap_open_score, gap_open_score, gap_open_score, gap_open_score, gap_open_score]);
    let gap_extend_8 = i16x8::from_array([gap_extend_score, gap_extend_score, gap_extend_score, gap_extend_score, gap_extend_score, gap_extend_score, gap_extend_score, gap_extend_score]);
    let left_mask_1 = i16x8::from_array([0, 1, 1, 1, 1, 1, 1, 1]);
    let right_mask_7 = i16x8::from_array([1, 0, 0, 0, 0, 0, 0, 0]);
    let mut HH: Vec<i16x8> = vec![i16x8::from_array([0, 0, 0, 0, 0, 0, 0, 0]); num_seq_vec];
    let mut EE: Vec<i16x8> = vec![i16x8::from_array([0, 0, 0, 0, 0, 0, 0, 0]); num_seq_vec];
    for i in 0..seq_x.len() {
        let mut X = i16x8::from_array([0, 0, 0, 0, 0, 0, 0, 0]);
        let mut F = i16x8::from_array([0, 0, 0, 0, 0, 0, 0, 0]);
        // get the database base
        let data_base_index = hash_table.get(&seq_x[i]).unwrap();
        for j in 0..num_seq_vec {
            let mut H = HH[j].clone();
            let mut E = EE[j].clone();
            //println!("H {:?}", H);
            let T1 = H.rotate_elements_left::<7>() * right_mask_7;
            //println!("T1 {:?}", T1);
            H = (H.rotate_elements_right::<1>() * left_mask_1) + X;
            //println!("H after rotate {:?}", H);
            X = T1;
            // convert MM to simd
            let MM_simd = MM_simd_full[*data_base_index][j];
            //let MM_simd = i16x8::from_array(MM[0..8].try_into().expect(""));
            //println!("MM {:?}", MM_simd);
            H = H + MM_simd;
            H = H.simd_max(E);

            F = F.rotate_elements_left::<7>() * right_mask_7;
            F = (H.rotate_elements_right::<1>() * left_mask_1) + F;
            // make simd of gap open gap extend 
            F = F - gap_extend_8 - gap_open_8;
            let mut T2 = F.clone();
            let all_non_positive = T2.simd_le(i16x8::splat(0)).all();
            if all_non_positive {
                F = H;
            }
            else {
                for _iter in 0..8 {
                    // lshift 2 t2, for gap extends
                    T2 = T2.rotate_elements_right::<1>() * left_mask_1;
                    F = F.simd_max(T2);
                }
                // update H based on max F
                H = H.simd_max(F);
                F = F + gap_open_8;
                F = F.simd_max(H);
            }
            
            // Final update
            HH[j] = H.clone();
            H = H - gap_open_8;
            E = E.simd_max(H);
            E = E - gap_extend_8;
            //println!("E at end {:?}", E);
            EE[j] = E;
            //println!("");
        }
        println!("{:?}", HH);
    }
}

pub fn fake_pairwise_simd (seq_x: &Vec<u8>, seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32) {
    println!("fake");
    // filling out score matrices and back matrix
    // the number of vectors for query SIMD 8 for now
    let num_seq_vec = seq_y.len() / 8;
    // make vectors of 8 of size num_seq_vec
    let mut HH = vec![vec![0; 8]; num_seq_vec];
    let mut EE = vec![vec![0; 8]; num_seq_vec];
    for i in 0..seq_x.len() {
        let mut X = 0;
        let mut F = vec![0; 8];
        for j in 0..num_seq_vec {
            let mut H = HH[j].clone();
            let mut E = EE[j].clone();
            //println!("H {:?}", H);
            let T1 = H[7];
            //println!("T1 {:?}", T1);
            H = vec![X, H[0], H[1], H[2], H[3], H[4], H[5], H[6]];
            //println!("H after rotate {:?}", H);
            X = T1;
            let mut MM = vec![0; 8];
            // make a match vector of 8, just for testing, for the real thing profiling is required
            let base_to_match = seq_x[i];
            for iter in 0..8 {
                if base_to_match == seq_y[(j * 8) + iter] {
                    MM[iter] = match_score;
                }
                else {
                    MM[iter] = mismatch_score;
                }
                // add match mismatch score to H
                H[iter] = H[iter] + MM[iter];
                // make H max of E or H
                if E[iter] > H[iter] {
                    H[iter] = E[iter];
                }
                // make F
            }
            //println!("MM {:?}", MM);
            F = vec![F[7], H[0], H[1], H[2], H[3], H[4], H[5], H[6]];
            for iter in 0..8 {
                F[iter] = F[iter] - gap_open_score - gap_extend_score;
            }
            let mut T2 = F.clone();
            for _iter in 0..8 {
                // lshift 2 t2, for gap extends
                T2 = vec![0, T2[0], T2[1], T2[2], T2[3], T2[4], T2[5], T2[6]];
                for iter2 in 0..8 {
                    // update if better
                    if T2[iter2] > F[iter2] {
                        F[iter2] = T2[iter2];
                    }
                }
            }
            // update H based on max F
            for iter in 0..8 {
                if F[iter] > H[iter] {
                    H[iter] = F[iter];
                }
                F[iter] = F[iter] + gap_open_score;
                if H[iter] > F[iter] {
                    F[iter] = H[iter];
                }
            }
            // Final update
            HH[j] = H.clone();
            for iter in 0..8 {
                if (H[iter] - gap_open_score) > E[iter] {
                    E[iter] = H[iter] - gap_open_score;
                }   
                E[iter] = E[iter] - gap_extend_score;
            }
            //println!("E at end {:?}", E);
            EE[j] = E;
            //println!("");
        }
        println!("{:?}", HH);
    }

}

pub fn pairwise (seq_x: &Vec<u8>, seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32, band_size: usize) -> (Vec<u8>, isize) {
    // variables to save results
    let mut align_vec: Vec<u8> = Vec::new();

    // make one matrix 
    let mut pair_wise_matrix: Vec<Vec<PairwiseMatrixCell>> = vec![vec![PairwiseMatrixCell { match_score: (0), del_score: (0), ins_score: (0), back: ('X') }; seq_y.len() + 1]; seq_x.len() + 1];

    // fill out the first row and colomn
    let temp_value = gap_open_score as isize + gap_extend_score as isize;
    //pair_wise_matrix[0][1] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('d') };
    //pair_wise_matrix[1][0] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('i') };
    for i in 2..seq_y.len() + 1 {
        //let temp_value = pair_wise_matrix[0][i - 1].del_score + gap_extend_score as isize;
        //pair_wise_matrix[0][i] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('d') };
    }
    for i in 2..seq_x.len() + 1 {
        //let temp_value = pair_wise_matrix[i - 1][0].ins_score + gap_extend_score as isize;
        //pair_wise_matrix[i][0] = PairwiseMatrixCell { match_score: (temp_value), del_score: (temp_value), ins_score: (temp_value), back: ('i') };
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