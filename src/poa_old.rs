// Copyright 2017-2025 Brett Bowman, Jeff Knaggs, Minindu Weerakoon
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of multiple homologous sequences.
//!
//! - time complexity: `O(N^2 * L^2)`, where `N` is the number of sequences and `L` is the length of each sequence.
//!
//! For the original concept and theory, see:
//! * Lee, Christopher, Catherine Grasso, and Mark F. Sharlow. "Multiple sequence alignment using
//! partial order graphs." Bioinformatics 18.3 (2002): 452-464.
//! * Lee, Christopher. "Generating consensus sequences from partial order multiple sequence
//! alignment graphs." Bioinformatics 19.8 (2003): 999-1008.
//!
//! For a modern reference implementation, see poapy:
//! https://github.com/ljdursi/poapy
//!
//! # Example
//!
//! ```
//! use bio::alignment::pairwise::Scoring;
//! use bio::alignment::poa::*;
//!
//! let x = b"AAAAAAA";
//! let y = b"AABBBAA";
//! let z = b"AABCBAA";
//!
//! let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
//! let mut aligner = Aligner::new(scoring, x);
//! // z differs from x in 3 locations
//! assert_eq!(aligner.global(z).alignment().score, 1);
//! aligner.global(y).add_to_graph();
//! // z differs from x and y's partial order alignment by 1 base
//! assert_eq!(aligner.global(z).alignment().score, 5);
//! ```

use std::cmp::{max};

use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs
pub type POAGraph = Graph<u8, i32, Directed, usize>;

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
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
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct Traceback {
    rows: usize,
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<(Vec<i32>, usize, usize)>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    fn with_capacity(m: usize, n: usize) -> Self {
        // each row of matrix contain start end position and vec of traceback cells
        let matrix: Vec<(Vec<i32>, usize, usize)> = vec![(vec![], 0, n + 1); m + 1];
        Traceback {
            rows: m,
            cols: n,
            last: NodeIndex::new(0),
            matrix,
        }
    }
    /// Populate the first row of the traceback matrix
    fn initialize_scores(&mut self, gap_open: i32, yclip: i32) {
        for j in 0..=self.cols {
            self.matrix[0].0.push(max((j as i32) * gap_open, yclip));
        }
        self.matrix[0].0[0] = 0;
    }

    fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            last: NodeIndex::new(0),
            matrix: Vec::new(),
        }
    }

    // create a new row according to the parameters
    fn new_row(
        &mut self,
        row: usize,
        size: usize,
        gap_open: i32,
        xclip: i32,
        start: usize,
        end: usize,
    ) {
        self.matrix[row].1 = start;
        self.matrix[row].2 = end;
        // when the row starts from the edge
        if start == 0 {
            self.matrix[row].0.push(max((row as i32) * gap_open, xclip));
        } else {
            self.matrix[row].0.push(MIN_SCORE);
        }
        for _ in 1..=size {
            self.matrix[row].0.push(MIN_SCORE);
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: i32) {
        // set the matrix cell if in band range
        if !(self.matrix[i].1 > j || self.matrix[i].2 < j) {
            let real_position = j - self.matrix[i].1;
            self.matrix[i].0[real_position] = cell;
        }
    }

    fn get(&self, i: usize, j: usize) -> i32 {
        // get the matrix cell if in band range else return the appropriate values
        if !(self.matrix[i].1 > j || self.matrix[i].2 <= j || self.matrix[i].0.is_empty()) {
            let real_position = j - self.matrix[i].1;
            self.matrix[i].0[real_position]
        }
        // out of band
        else {
            MIN_SCORE
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
#[derive(Default, Clone, Debug)]
pub struct Aligner {
    pub traceback: Traceback,
    query: Vec<u8>,
    pub poa: Poa,
}

impl  Aligner {
    /// Create new instance.
    pub fn new(match_score: i32, mismatch_score: i32, gap_open_score: i32, reference: &Vec<u8>, x_clip: i32, y_clip: i32, _band_size: i32) -> Self {
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(match_score, mismatch_score, gap_open_score, x_clip, y_clip, reference),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self) -> (Vec<u8>, Vec<usize>) {
        let alignment = self.poa.recalculate_alignment(&self.traceback);
        println!("\n score old aligner {}", alignment.score);
        let path_indices = self.poa.add_alignment(&alignment, &self.query);
        let mut path_bases = vec![];
        for path_index in &path_indices {
            path_bases.push(self.poa.graph.raw_nodes()[*path_index].weight);
        }
        (path_bases, path_indices)
    }

    /// Custom align a given query against the graph with custom xclip and yclip penalties.
    pub fn custom(&mut self, query: &Vec<u8>) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);
        self
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

/// A partially ordered alignment graph
///
/// A directed acyclic graph datastructure that represents the topology of a
/// traceback matrix.
#[derive(Default, Clone, Debug)]
pub struct Poa {
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    x_clip: i32,
    y_clip: i32,
    pub graph: POAGraph,
}

impl Poa {
    /// Create a new POA graph from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    pub fn from_string(match_score: i32, mismatch_score: i32, gap_open_score: i32, x_clip: i32, y_clip: i32, seq: &Vec<u8>) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        Poa { match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score, x_clip: x_clip, y_clip: y_clip, graph}
    }
    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    pub fn custom(&self, query: &Vec<u8>) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        // save score location of the max scoring node for the query for suffix clipping
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.gap_open_score, self.y_clip);
        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1; // 0 index is for initialization so we start at 1
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            traceback.new_row(
                i,
                n + 1,
                self.gap_open_score,
                self.x_clip,
                0,
                n + 1,
            );
            // query base and its index in the DAG (traceback matrix rows)
            for (query_index, query_base) in query.iter().enumerate() {
                let j = query_index + 1; // 0 index is initialized so we start at 1
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    if r == *query_base {
                        traceback.get(0, j - 1) + self.match_score
                    }
                    else {
                        traceback.get(0, j - 1) + self.mismatch_score
                    }
                } else {
                    let mut max_cell = max(MIN_SCORE, self.x_clip);
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        let temp_score;
                        if r == *query_base {
                            temp_score = self.match_score;
                        }
                        else {
                            temp_score = self.mismatch_score;
                        }
                        max_cell = max(
                            max_cell,
                            max(traceback.get(i_p, j - 1) + temp_score, traceback.get(i_p, j) + self.gap_open_score)
                        );
                    }
                    max_cell
                };
                let score = max(max_cell, traceback.get(i, j - 1)+ self.gap_open_score);
                traceback.set(i, j, score);
            }
        }
        traceback
    }
    // Recalculate the alignment using the traceback
    // We need this in poa because we have to use the scoring 
    pub fn recalculate_alignment(&self, traceback: &Traceback) -> Alignment {
        // Get the alignment by backtracking and recalculating stuff
        let mut ops: Vec<AlignmentOperation> = vec![];
        // loop until we reach a node with no incoming
        let mut curr_node = traceback.last.index() + 1;
        let mut curr_query = traceback.cols;
        let final_score = traceback.get(curr_node, curr_query);
        println!("old alignment score {}", final_score);
        loop {
            let mut current_alignment_operation = AlignmentOperation::Match(None);
            let current_cell_score = traceback.get(curr_node, curr_query);
            let mut next_jump = 0;
            let mut next_node = 1;
            // check left if gap open difference with left
            let prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(NodeIndex::new(curr_node - 1), Incoming).collect();
            if current_cell_score == traceback.get(curr_node, curr_query - 1) + self.gap_open_score {
                current_alignment_operation = AlignmentOperation::Ins(Some(curr_node - 1));
                next_jump = curr_query - 1;
                next_node = curr_node;
            }
            else {
                for prev in &prevs {
                    let prev_node = prev.index() + 1;
                    // Top
                    if current_cell_score == traceback.get(prev_node, curr_query) + self.gap_open_score {
                        current_alignment_operation = AlignmentOperation::Del(Some((0, curr_node - 1)));
                        next_jump = curr_query;
                        next_node = prev_node;
                    }
                    // Diagonal
                    else if (current_cell_score == traceback.get(prev_node, curr_query - 1) + self.match_score) ||
                    (current_cell_score == traceback.get(prev_node, curr_query - 1) + self.mismatch_score) {
                        current_alignment_operation = AlignmentOperation::Match(Some((prev_node - 1, curr_node - 1)));
                        next_jump = curr_query - 1;
                        next_node = prev_node;
                    }
                }
            }
            // if out of band we can go left and if hit an edge go up
            ops.push(current_alignment_operation);
            // iterate to next
            curr_query = next_jump;
            curr_node = next_node;
            // break point
            if prevs.len() == 0 || curr_query == 0 {
                //if at end but not at start of query add bunch of ins(None)
                if prevs.len() == 0 {
                    if curr_query > 0 {
                        for _ in 0..curr_query {
                            ops.push(AlignmentOperation::Ins(None));
                        }
                    }
                } else {
                    // push del until we hit no prevs
                    loop {
                        let prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(NodeIndex::new(curr_node - 1), Incoming).collect();
                        if prevs.len() == 0 {break}
                        ops.push(AlignmentOperation::Del(Some((0, prevs[0].index()))));
                        curr_node = prevs[0].index() + 1;
                    }
                }
                break;
            }
        }
        ops.reverse();
        Alignment {
            score: final_score as i32,
            operations: ops
        }
    }
    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    pub fn add_alignment(&mut self, aln: &Alignment, seq: &Vec<u8>) -> Vec<usize> {
        let mut path_indices = vec![];
        let head = Topo::new(&self.graph).next(&self.graph).unwrap();
        let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
        let mut i: usize = 0;
        let mut edge_not_connected: bool = false;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    let node: NodeIndex<usize> = NodeIndex::new(head.index());
                    if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        path_indices.push(node.index());
                        if edge_not_connected {
                            self.graph.add_edge(prev, node, 1);
                        }
                        edge_not_connected = false;
                        prev = node;
                    }
                    else {
                        path_indices.push(node.index());
                    }
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        edge_not_connected = false;
                    }
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        path_indices.push(node.index());
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        path_indices.push(node.index());
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
                    path_indices.push(node.index());
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
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                // we should only have to skip over deleted nodes and xclip
                AlignmentOperation::Del(Some((_, x))) => {
                    path_indices.push(*x);
                } 
                AlignmentOperation::Del(None) => {} 
                AlignmentOperation::Xclip(_) => {}
                AlignmentOperation::Yclip(_, r) => {
                    i = *r;
                }
            }
        }
        path_indices
    }
}