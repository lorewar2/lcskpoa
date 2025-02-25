use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;
use petgraph::{Directed, Graph, Incoming};
pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs
pub type POAGraph = Graph<u8, i32, Directed, usize>;
use std::path;
use std::simd::i32x8;
use std::simd::cmp::SimdOrd;
use std::collections::HashMap;

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
    Xclip(usize),
    Yclip(usize, usize), // to, from
}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Alignment {
    pub score: i32,
    operations: Vec<AlignmentOperation>,
}

#[derive(Default, Clone, Debug)]
pub struct  SimdTracker {
    last_node: usize,
    last_query: usize,
    gap_open: i32, // required for sending fake data of unbanded sections
    simd_matrix: Vec<Vec<i32x8>>, // each row has the start end info and simd vecs (start end info is for simd vec indices)
    start_end_tracker: Vec<(usize, usize)>,
}

impl SimdTracker {
    fn new (num_of_nodes: usize, query_len: usize, simd_vecs: usize, gap_open: i32) -> Self {
        let mut simd_matrix: Vec<Vec<i32x8>> = vec![vec![]; num_of_nodes];
        let start_end_tracker = vec![(0, simd_vecs); num_of_nodes];
        let gap_open_8 = i32x8::splat(-gap_open);
        // make the index 0 one
        let gap_multiplier = i32x8::from_array([1, 2, 3, 4, 5, 6, 7, 8]);
        simd_matrix[0] = (0..simd_vecs).map(|j| {
            let base_offset = (j * 8) as i32;
            (gap_multiplier + i32x8::splat(base_offset)) * -gap_open_8
        }).collect();
        SimdTracker {
            last_node: num_of_nodes,
            last_query: query_len,
            gap_open: gap_open,
            simd_matrix: simd_matrix,
            start_end_tracker: start_end_tracker,
        }
    }

    fn with_capacity(num_of_nodes: usize, query_len: usize, simd_vecs: usize, gap_open: i32) -> Self {
        let start_end_tracker = vec![(0, simd_vecs); num_of_nodes];
        let mut simd_matrix: Vec<Vec<i32x8>> = Vec::with_capacity(num_of_nodes);
        let gap_open_8 = i32x8::splat(-gap_open);
        // make the index 0 one
        let gap_multiplier = i32x8::from_array([1, 2, 3, 4, 5, 6, 7, 8]);
        let min_score_8 = i32x8::splat(MIN_SCORE);
        // initialize first row
        simd_matrix.push((0..simd_vecs).map(|j| {
            let base_offset = (j * 8) as i32;
            (gap_multiplier + i32x8::splat(base_offset)) * -gap_open_8
        }).collect());
        for _i in 1..num_of_nodes {
            simd_matrix.push(vec![min_score_8; simd_vecs]);
        }
        // make the other rows as well to test
        SimdTracker {
            last_node: num_of_nodes,
            last_query: query_len,
            gap_open: gap_open,
            simd_matrix: simd_matrix,
            start_end_tracker: start_end_tracker,
        }
    }
    // allocate the matrix row with MIN SCORE stuff
    fn new_row(&mut self, row: usize, start: usize, end: usize){
        self.start_end_tracker[row] = (start, end);
        for _ in start..end + 1{
            self.simd_matrix[row].push(i32x8::splat(MIN_SCORE));
        }
        
    }
    // get function, if not in band do something, try to get it back to band
    fn get(&self, i: usize, j: usize) -> i32x8 {
        // get the matrix cell if in band range else return the appropriate values
        if !(self.start_end_tracker[i].0 > j || self.start_end_tracker[i].1 <= j || self.simd_matrix[i].is_empty()) {
            let real_position = j - self.start_end_tracker[i].0;
            self.simd_matrix[i][real_position]
        }
        // this should happen, but if it did try to control
        else if j > self.start_end_tracker[i].1 {
            let neg_10 = i32x8::splat(-(10 * (1_000_000 - i as i32))); //will not go up or diagonal, **gap end/mismatch should not be 10** test this :C
            let gap_open = i32x8::splat(self.gap_open as i32);
            let j_multi = i32x8::splat(j as i32 * 8);
            let j_increment = i32x8::from_array([1, 2, 3, 4, 5, 6, 7, 8]);
            //println!("Go left");
            return ((j_multi + j_increment) * gap_open) + neg_10;
        }
        // make it go up to meet the band, modify the values do it does not go left should be ok with same numbers
        else {
            let neg_10 = i32x8::splat(self.gap_open * (1_000_000 - i as i32));
            //println!("Go UP");
            return i32x8::from_array([1, 1, 1, 1, 1, 1, 1, 1]) * neg_10;
        }
    }
    fn set(&mut self, i: usize, j: usize, simd: i32x8) {
        // set the matrix cell if in band range
        if !(self.start_end_tracker[i].0 > j || self.start_end_tracker[i].1 < j) {
            let real_position = j - self.start_end_tracker[i].0;
            self.simd_matrix[i][real_position] = simd;
        }
    }
}

pub struct Aligner {
    query: Vec<u8>,
    poa: Poa
}

impl Aligner {
    /// Create new instance.
    pub fn new(match_score: i32, mismatch_score: i32, gap_open_score: i32, reference: &Vec<u8>) -> Self {
        Aligner {
            query: reference.to_vec(),
            poa: Poa::from_string(match_score, mismatch_score, gap_open_score, reference),
        }
    }

    pub fn global_simd(&mut self, query: &Vec<u8>) -> (Vec<u8>, Vec<usize>) {
        self.query = query.to_vec();
        let simd_tracker = self.poa.custom_simd(query);
        let alignment = self.poa.recalculate_alignment(simd_tracker);
        let path_indices = self.poa.add_alignment(&alignment, &self.query);
        let mut path_bases = vec![];
        for path_index in &path_indices {
            path_bases.push(self.poa.graph.raw_nodes()[*path_index].weight);
        }
        (path_bases, path_indices)
    }

    pub fn global_simd_banded(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, band_size: usize) -> (Vec<u8>, Vec<usize>) {
        self.query = query.to_vec();
        let simd_tracker = self.poa.custom_simd_banded(query, lcsk_path, band_size);
        let alignment = self.poa.recalculate_alignment(simd_tracker);
        let path_indices = self.poa.add_alignment(&alignment, &self.query);
        let mut path_bases = vec![];
        for path_index in &path_indices {
            path_bases.push(self.poa.graph.raw_nodes()[*path_index].weight);
        }
        (path_bases, path_indices)
    }

    /// Return alignment graph.
    pub fn graph(&self) -> &POAGraph {
        &self.poa.graph
    }

    /// Return the consensus sequence generated from the POA graph.
    pub fn consensus(&self) -> Vec<u8> {
        let mut consensus: Vec<u8> = vec![];
        let max_index = self.poa.graph.node_count();
        let mut weight_score_next_vec: Vec<(i32, i32, usize)> = vec![(0, 0, 0); max_index + 1];
        let mut topo = Topo::new(&self.poa.graph);
        // go through the nodes topologically
        while let Some(node) = topo.next(&self.poa.graph) {
            let mut best_weight_score_next: (i32, i32, usize) = (0, 0, usize::MAX);
            let neighbour_nodes = self.poa.graph.neighbors_directed(node, Incoming);
            // go through the incoming neighbour nodes
            for neighbour_node in neighbour_nodes {
                let neighbour_index = neighbour_node.index();
                let neighbour_score = weight_score_next_vec[neighbour_index].1;
                let edges = self.poa.graph.edges_connecting(neighbour_node, node);
                let weight = edges.map(|edge| edge.weight()).sum();
                let current_node_score = weight + neighbour_score;
                // save the neighbour node with the highest weight and score as best
                if (weight, current_node_score, neighbour_index) > best_weight_score_next {
                    best_weight_score_next = (weight, current_node_score, neighbour_index);
                }
            }
            weight_score_next_vec[node.index()] = best_weight_score_next;
        }
        // get the index of the max scored node (end of consensus)
        let mut pos = weight_score_next_vec
            .iter()
            .enumerate()
            .max_by_key(|(_, &value)| value.1)
            .map(|(idx, _)| idx)
            .unwrap();
        // go through weight_score_next_vec appending to the consensus
        while pos != usize::MAX {
            consensus.push(self.poa.graph.raw_nodes()[pos].weight);
            pos = weight_score_next_vec[pos].2;
        }
        consensus.reverse();
        consensus
    }
}

#[derive(Default, Clone, Debug)]
pub struct Poa {
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    pub graph: POAGraph,
    pub memory_usage: usize,
}

impl Poa {
    pub fn from_string(match_score: i32, mismatch_score: i32, gap_open_score: i32, seq: &Vec<u8>) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        Poa { match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score, graph, memory_usage: 0}
    }
    pub fn profile_query (seq_y: &Vec<u8>, match_score: i32, mismatch_score: i32) -> Vec<Vec<i32x8>> {
        let num_seq_vec = (seq_y.len() + 7) / 8;
        let mut mm_simd = vec![vec![]; 4];
        // make 4 vectors for query
        let score_table = [0, mismatch_score];
        // go through the query and populate the entries
        for index_simd in 0..num_seq_vec {
            let mut temp_a = [0; 8];
            let mut temp_c = [0; 8];
            let mut temp_g = [0; 8];
            let mut temp_t = [0; 8];
            let start = index_simd * 8;
            let end = (start + 8).min(seq_y.len());
            let padded_seq = seq_y[start..end]
                .iter()
                .chain(std::iter::repeat(&0))
                .take(8);
            for (i, &base) in padded_seq.enumerate() {
                temp_a[i] = if base == 65 { match_score } else { score_table[(base > 0) as usize] };
                temp_c[i] = if base == 67 { match_score } else { score_table[(base > 0) as usize] };
                temp_g[i] = if base == 71 { match_score } else { score_table[(base > 0) as usize] };
                temp_t[i] = if base == 84 { match_score } else { score_table[(base > 0) as usize] };
            }
            mm_simd[0].push(i32x8::from_array(temp_a));
            mm_simd[1].push(i32x8::from_array(temp_c));
            mm_simd[2].push(i32x8::from_array(temp_g));
            mm_simd[3].push(i32x8::from_array(temp_t));
        }
        mm_simd
    }
    pub fn custom_simd_banded (&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, band_size: usize) -> SimdTracker {
        let mut hash_table = HashMap::new();
        hash_table.insert(65, 0);
        hash_table.insert(67, 1);
        hash_table.insert(71, 2);
        hash_table.insert(84, 3);
        let mm_simd_full = Poa::profile_query(query, self.match_score, self.mismatch_score);
        // other simd stuff required
        let gap_open_score = -self.gap_open_score as i32;
        let gap_open_8 = i32x8::splat(gap_open_score);
        let zero_8 = i32x8::splat(0);
        let min_score_8 = i32x8::splat(MIN_SCORE);
        let left_mask_1 = i32x8::from_array([0, 1, 1, 1, 1, 1, 1, 1]);
        let right_mask_7 = i32x8::from_array([1, 0, 0, 0, 0, 0, 0, 0]);
        assert!(self.graph.node_count() != 0);
        let query_len = query.len();
        // dimensions of the traceback matrix
        let num_seq_vec = (query_len as f64 / 8.0).ceil() as usize;
        // use lcsk path here for making matrix
        let mut simd_tracker = SimdTracker::new(self.graph.node_count(), query_len, num_seq_vec, self.gap_open_score);
        // construct the score matrix (O(n^2) space)
        let mut index = 0;
        let mut topo = Topo::new(&self.graph);
        // THIS STUFF FOR BAND PART 1 
        let mut no_kmers = false;
        if lcsk_path.len() == 0 {
            no_kmers = true;
        }
        println!("no kmers {}", no_kmers);
        let mut start_banding_query_node = (0, 0);
        let mut end_banding_query_node = &(0, 0);
        let mut banding_started = false;
        let mut banding_ended = false;
        let mut current_lcsk_path_index = 0;
        if !no_kmers {
            start_banding_query_node = lcsk_path[0];
            end_banding_query_node = lcsk_path.last().unwrap();
        }
        while let Some(node) = topo.next(&self.graph) {
            let mut f = zero_8;
            f[7] = (index + 1) * -gap_open_score;
            let r = self.graph.raw_nodes()[node.index()].weight;
            let i = node.index(); 
            let mut start = 0;
            let mut end = query_len;
            if !no_kmers {
                if banding_started == false {
                    end = start_banding_query_node.0 + band_size;
                }
                else if banding_ended == true {
                    start = if band_size > end_banding_query_node.0 {
                        0
                    } else {
                        end_banding_query_node.0 - band_size
                    };
                }
                else{
                    start = if band_size > lcsk_path[current_lcsk_path_index].0 {
                        0
                    } else {
                        lcsk_path[current_lcsk_path_index].0 - band_size
                    };
                    if lcsk_path.len() < current_lcsk_path_index + 1 {
                        end = lcsk_path[current_lcsk_path_index + 1].0 + band_size;
                    }
                    else {
                        end = lcsk_path[current_lcsk_path_index].0 + band_size;
                    }
                }
                if banding_ended != true {
                    if lcsk_path[current_lcsk_path_index].1 == i {
                        current_lcsk_path_index += 1;
                    }
                    if start_banding_query_node.1 == i {
                        banding_started = true;
                    }
                    if end_banding_query_node.1 == i {
                        banding_ended = true;
                    }
                }
            }
            if end > query_len {
                end = query_len;
            }
            let start_simd = start / 8;
            let end_simd = (end / 8) + 1;
            if i != 0 {
                simd_tracker.new_row(i, start_simd, end_simd); 
            }
            simd_tracker.last_node = i;
            let data_base_index = hash_table.get(&r).unwrap();
            let mut prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(node, Incoming).collect();
            if prevs.len() == 0 {
                prevs.push(node);
            }
            // vertical and diagonal
            for prev_node in &prevs {
                let i_p: usize = prev_node.index(); // index of previous node
                let mut x = zero_8;
                x[0] = (index) * -gap_open_score;
                for simd_index in 0..num_seq_vec {
                    if start_simd > simd_index {
                        continue;
                    }
                    if end_simd < simd_index {
                        break;
                    }
                    let h_prev = simd_tracker.get(i_p, simd_index);
                    let mut h_curr;
                    if i_p == i {
                        h_curr = min_score_8; 
                    }
                    else {
                        h_curr = simd_tracker.get(i, simd_index);
                    }
                    let e = h_prev - gap_open_8;
                    let mm_simd = mm_simd_full[*data_base_index][simd_index];
                    let t1 = h_prev.rotate_elements_left::<7>() * right_mask_7;
                    let mut t2 = (h_prev.rotate_elements_right::<1>() * left_mask_1) + x;
                    x = t1;
                    t2 = t2 + mm_simd;
                    h_curr = h_curr.simd_max(t2);
                    h_curr = h_curr.simd_max(e);
                    simd_tracker.set(i, simd_index, h_curr);
                }
            }
            // horizontal NEeds fixing, non simd is faster here
            for simd_index in 0..num_seq_vec {
                if start_simd > simd_index {
                    continue;
                }
                if end_simd < simd_index {
                    break;
                }
                let h = simd_tracker.get(i, simd_index);
                let mut t3 = f[7];
                let mut max_vec = zero_8;
                for iter in 0..8 {
                    let temp = h[iter];
                    t3 = t3 - gap_open_score;
                    if temp > t3 {
                        t3 = temp;
                    }
                    max_vec[iter] = t3;
                }
                simd_tracker.set(i, simd_index, max_vec);
                f = max_vec;
            }
            index += 1;
        }
        simd_tracker
    }

    pub fn custom_simd(&mut self, query: &Vec<u8>) -> SimdTracker {
        let mut hash_table = HashMap::new();
        hash_table.insert(65, 0);
        hash_table.insert(67, 1);
        hash_table.insert(71, 2);
        hash_table.insert(84, 3);
        let mm_simd_full = Poa::profile_query(query, self.match_score, self.mismatch_score);
        // other simd stuff required
        let gap_open_score = -self.gap_open_score as i32;
        let gap_open_8 = i32x8::splat(gap_open_score);
        let zero_8 = i32x8::splat(0);
        let min_score_8 = i32x8::splat(MIN_SCORE);
        let left_mask_1 = i32x8::from_array([0, 1, 1, 1, 1, 1, 1, 1]);
        let right_mask_7 = i32x8::from_array([1, 0, 0, 0, 0, 0, 0, 0]);
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let query_len  =  query.len();
        let num_seq_vec = (query_len as f64 / 8.0).ceil() as usize;
        let mut simd_tracker = SimdTracker::new(self.graph.node_count(), query_len, num_seq_vec, self.gap_open_score);
        // construct the score matrix (O(n^2) space)
        let mut index = 0;
        let mut topo = Topo::new(&self.graph);
        // required stuff for backtrace
        while let Some(node) = topo.next(&self.graph) {
            let mut f = zero_8;
            f[7] = (index + 1) * -gap_open_score;
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight;
            let i = node.index(); 
            let start_simd = 0;
            let end_simd = (query_len / 8) + 1;
            if i != 0 {
                simd_tracker.new_row(i, start_simd, end_simd); 
            }
            simd_tracker.last_node = i;
            let data_base_index = hash_table.get(&r).unwrap();
            let mut prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(node, Incoming).collect();
            if prevs.len() == 0 {
                prevs.push(node);
            }
            // vertical and diagonal
            for prev_node in &prevs {
                let i_p: usize = prev_node.index(); // index of previous node
                let mut x = zero_8;
                x[0] = (index) * -gap_open_score;
                for simd_index in 0..num_seq_vec {
                    let h_prev = simd_tracker.get(i_p, simd_index);
                    let mut h_curr;
                    if i_p == i {
                        h_curr = min_score_8; 
                    }
                    else {
                        h_curr = simd_tracker.get(i, simd_index);
                    }
                    let e = h_prev - gap_open_8;
                    let mm_simd = mm_simd_full[*data_base_index][simd_index];
                    let t1 = h_prev.rotate_elements_left::<7>() * right_mask_7;
                    let mut t2 = (h_prev.rotate_elements_right::<1>() * left_mask_1) + x;
                    x = t1;
                    t2 = t2 + mm_simd;
                    h_curr = h_curr.simd_max(t2);
                    h_curr = h_curr.simd_max(e);
                    simd_tracker.set(i, simd_index, h_curr);
                }
            }
            // horizontal NEeds fixing, non simd is faster here
            for simd_index in 0..num_seq_vec {
                let h = simd_tracker.get(i, simd_index);
                let mut t3 = f[7];
                let mut max_vec = zero_8;
                for iter in 0..8 {
                    let temp = h[iter];
                    t3 = t3 - gap_open_score;
                    if temp > t3 {
                        t3 = temp;
                    }
                    max_vec[iter] = t3;
                }
                simd_tracker.set(i, simd_index, max_vec);
                f = max_vec;
                //print!("{:?} ", max_vec);
            }
            //println!("");
            index += 1;
        }
        simd_tracker
    }

    pub fn recalculate_alignment (&mut self, simd_tracker: SimdTracker) -> Alignment {
        // Get the alignment by backtracking and recalculating stuff
        let mut ops: Vec<AlignmentOperation> = vec![];
        // loop until we reach a node with no incoming
        let mut current_node = simd_tracker.last_node;
        let mut current_query = simd_tracker.last_query - 1;
        
        let simd_index = (current_query) / 8;
        let simd_inner_index = (current_query) % 8;
        let simd_vec_obtained = simd_tracker.get(current_node, simd_index);
        let final_score = simd_vec_obtained[simd_inner_index];
        println!("simd score: {} vec {:?} {} {} {:?}", final_score, simd_vec_obtained, simd_index, simd_inner_index, simd_tracker.start_end_tracker[current_node]);
        loop {
            let mut current_alignment_operation = AlignmentOperation::Match(None);
            //check the score ins left of query
            let prev_simd_index = (current_query - 1) / 8;
            let prev_simd_inner_index = (current_query - 1) % 8;

            let simd_index = (current_query) / 8;
            let simd_inner_index = (current_query) % 8;

            let simd_vec_obtained = simd_tracker.get(current_node, simd_index);
            let current_cell_score = simd_vec_obtained[simd_inner_index];
            let mut next_jump = 0;
            let mut next_node = 1;
            // Check left if gap open difference with left
            let prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(NodeIndex::new(current_node), Incoming).collect();
            let simd_prev_vec_obtained = simd_tracker.get(current_node, prev_simd_index);
            if current_cell_score == simd_prev_vec_obtained[prev_simd_inner_index] + self.gap_open_score {
                current_alignment_operation = AlignmentOperation::Ins(Some(current_node));
                next_jump = current_query - 1;
                next_node = current_node;
            }
            else {
                for prev in &prevs {
                    let i_p = prev.index();
                    // Top
                    let simd_vec_obtained = simd_tracker.get(i_p, simd_index);
                    let simd_prev_vec_obtained = simd_tracker.get(i_p, prev_simd_index);
                    if current_cell_score == simd_vec_obtained[simd_inner_index] + self.gap_open_score {
                        //current_alignment_operation = AlignmentOperation::Del(None);
                        current_alignment_operation = AlignmentOperation::Del(Some((0, current_node)));
                        next_jump = current_query;
                        next_node = i_p;
                    }
                    // Diagonal
                    else if (current_cell_score == simd_prev_vec_obtained[prev_simd_inner_index] + self.match_score as i32) || (current_cell_score == simd_prev_vec_obtained[prev_simd_inner_index] + self.mismatch_score as i32) {
                        current_alignment_operation = AlignmentOperation::Match(Some((i_p, current_node)));
                        next_jump = current_query - 1;
                        next_node = i_p;
                    }
                }
            }
            ops.push(current_alignment_operation);
            current_query = next_jump;
            current_node = next_node;
            // Break point
            if prevs.len() == 0 || current_query == 0 {
                //if at end but not at start of query add bunch of ins(None)
                if prevs.len() == 0 {
                    if current_query > 0 {
                        for _ in 0..current_query {
                            ops.push(AlignmentOperation::Ins(None));
                        }
                    }
                } else {
                    // push del until we hit no prevs
                    loop {
                        let prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(NodeIndex::new(current_query - 1), Incoming).collect();
                        if prevs.len() == 0 {break}
                        ops.push(AlignmentOperation::Del(Some((0, prevs[0].index()))));
                        current_query = prevs[0].index();
                    }
                }
                ops.push(AlignmentOperation::Match(None));
                break;
            }
        }
        ops.reverse();
        Alignment {
            score: final_score as i32,
            operations: ops
        }
    }

    pub fn add_alignment(&mut self, aln: &Alignment, seq: &Vec<u8>) -> Vec<usize> {
        let mut path_indices = vec![];
        let head = Topo::new(&self.graph).next(&self.graph).unwrap();
        let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
        let mut i: usize = 0;
        let mut edge_not_connected: bool = false;
        let mut prev_is_match_none = false;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    let node: NodeIndex<usize> = NodeIndex::new(head.index());
                    if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        if !prev_is_match_none {
                            path_indices.push(node.index());
                        }   
                        println!("match none 1 {}", node.index());
                        if edge_not_connected {
                            self.graph.add_edge(prev, node, 1);
                        }
                        edge_not_connected = false;
                        prev = node;
                    }
                    else {
                        if !prev_is_match_none {
                            path_indices.push(node.index());
                        }   
                        println!("match none 2 {}", node.index());
                    }
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        edge_not_connected = false;
                    }
                    i += 1;
                    prev_is_match_none = true;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        path_indices.push(node.index());
                        println!("match some 1 {}", node.index());
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        path_indices.push(node.index());
                        println!("match some 2 {}", node.index());
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                if prev.index() != head.index() && prev.index() != node.index() {
                                    self.graph.add_edge(prev, node, 1);
                                }
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    let node = self.graph.add_node(seq[i]);
                    if !prev_is_match_none {
                        path_indices.push(node.index());
                    } 
                    println!("Ins None 1 {}", node.index());
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                    }
                    prev = node;
                    edge_not_connected = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    path_indices.push(node.index());
                    println!("Ins Some 1 {}", node.index());
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                // we should only have to skip over deleted nodes and xclip
                AlignmentOperation::Del(Some((_, x))) => {
                    path_indices.push(*x);
                    println!("Del Some 1 {}", x);
                } 
                AlignmentOperation::Del(None) => {} 
                AlignmentOperation::Xclip(_) => {}
                AlignmentOperation::Yclip(_, r) => {
                    i = *r;
                }
            }
        }
        println!("{:?}", path_indices);
        path_indices
    }
}