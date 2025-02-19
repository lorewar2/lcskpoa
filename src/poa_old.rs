// Copyright 2017-2024 Brett Bowman, Jeff Knaggs, Minindu Weerakoon
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
use std::cmp::{max, Ordering};
use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;
use std::mem;
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

#[derive(Copy, Clone, Debug)]
pub struct TracebackCell {
    score: i32,
    op: AlignmentOperation,
}

impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.score == other.score
    }
}

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct Traceback {
    rows: usize,
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<(Vec<TracebackCell>, usize, usize)>,
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
        let matrix: Vec<(Vec<TracebackCell>, usize, usize)> = vec![(vec![], 0, n + 1); m + 1];
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
            self.matrix[0].0.push(max(
                TracebackCell {
                    score: (j as i32) * gap_open,
                    op: AlignmentOperation::Ins(None),
                },
                TracebackCell {
                    score: yclip,
                    op: AlignmentOperation::Yclip(0, j),
                },
            ));
        }
        self.matrix[0].0[0] = TracebackCell {
            score: 0,
            op: AlignmentOperation::Match(None),
        };
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
            self.matrix[row].0.push(max(
                TracebackCell {
                    score: (row as i32) * gap_open,
                    op: AlignmentOperation::Del(None),
                },
                TracebackCell {
                    score: xclip,
                    op: AlignmentOperation::Xclip(0),
                },
            ));
        } else {
            self.matrix[row].0.push(TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            });
        }
        for _ in 1..=size {
            self.matrix[row].0.push(TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            });
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: TracebackCell) {
        // set the matrix cell if in band range
        if !(self.matrix[i].1 > j || self.matrix[i].2 < j) {
            let real_position = j - self.matrix[i].1;
            self.matrix[i].0[real_position] = cell;
        }
    }

    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        // get the matrix cell if in band range else return the appropriate values
        if !(self.matrix[i].1 > j || self.matrix[i].2 <= j) && (self.matrix[i].0.len() > 0) {
            //println!("IN BAND, OK");
            let real_position = j - self.matrix[i].1;
            return &self.matrix[i].0[real_position];
        }
        // behind the band, met the edge
        else if j == 0 {
            //println!("BEHIND OF BAND, EDGE!!!");
            return &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Del(None),
            };
        }
        // infront of the band
        else if j >= self.matrix[i].2 {
            //println!("INFRONT OF BAND!!!");
            return &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Ins(None),
            };
        }
        // behind the band
        else {
            //println!("BEHIND BAND!!!");
            return &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            };
        }
    }

    pub fn alignment(&self) -> Alignment {
        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        // Now backtrack through the matrix to construct an optimal path
        let mut i = self.last.index() + 1;
        let mut j = self.cols;

        while i > 0 || j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.get(i, j).op.clone());
            match self.get(i, j).op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Match(None) => {
                    i = 0;
                    j -= 1;
                }
                AlignmentOperation::Del(None) => {
                    i -= 1;
                }
                AlignmentOperation::Ins(None) => {
                    j -= 1;
                }
                AlignmentOperation::Xclip(r) => {
                    i = r;
                }
                AlignmentOperation::Yclip(r, _) => {
                    j = r;
                }
            }
        }

        ops.reverse();
        //println!("old one {:?}", ops);
        Alignment {
            score: self.get(self.last.index() + 1, self.cols).score,
            operations: ops,
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
#[derive(Default, Clone, Debug)]
pub struct Aligner {
    traceback: Traceback,
    query: Vec<u8>,
    pub poa: Poa,
}

impl Aligner{
    /// Create new instance.
    pub fn empty (match_score: i32, mismatch_score: i32, gap_open_score: i32, x_clip: i32, y_clip: i32, _band_size: i32) -> Self {
        let reference = &vec![65];
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(match_score, mismatch_score, gap_open_score, x_clip, y_clip, reference),
        }
    }
    pub fn new(match_score: i32, mismatch_score: i32, gap_open_score: i32, reference: &Vec<u8>, x_clip: i32, y_clip: i32, _band_size: &i32) -> Self {
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(match_score, mismatch_score, gap_open_score, x_clip, y_clip, reference),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self) -> &mut Self {
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, &self.query);
        self
    }

    /// Return alignment of last added query against the graph.
    pub fn alignment(&self) -> Alignment {
        self.traceback.alignment()
    }

    /// Globally align a given query against the graph.
    pub fn global(&mut self, query: &Vec<u8>) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.x_clip,
            self.poa.y_clip,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.x_clip = MIN_SCORE;
        self.poa.y_clip = MIN_SCORE;

        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);

        // Set the clip penalties to the original values
        self.poa.x_clip = clip_penalties[0];
        self.poa.y_clip = clip_penalties[1];

        self
    }

    /// Semi-globally align a given query against the graph.
    pub fn semiglobal(&mut self, query: &Vec<u8>) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.x_clip,
            self.poa.y_clip,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.x_clip = MIN_SCORE;
        self.poa.y_clip = 0;

        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);

        // Set the clip penalties to the original values
        self.poa.x_clip = clip_penalties[0];
        self.poa.y_clip = clip_penalties[1];
        self
    }

    /// Locally align a given query against the graph.
    pub fn local(&mut self, query: &Vec<u8>) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.x_clip,
            self.poa.y_clip,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.x_clip = 0;
        self.poa.y_clip = 0;

        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);

        // Set the clip penalties to the original values
        self.poa.x_clip = clip_penalties[0];
        self.poa.y_clip = clip_penalties[1];

        self
    }

    /// Custom align a given query against the graph with custom xclip and yclip penalties.
    pub fn custom(&mut self, query: &Vec<u8>) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.custom(query);
        self
    }
    /// Globally align a given query against the graph.
    pub fn global_banded(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, bandwidth: usize) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.x_clip,
            self.poa.y_clip,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.x_clip = MIN_SCORE;
        self.poa.y_clip = MIN_SCORE;

        self.query = query.to_vec();
        self.traceback = self.poa.custom_banded(query, lcsk_path, bandwidth);

        // Set the clip penalties to the original values
        self.poa.x_clip = clip_penalties[0];
        self.poa.y_clip = clip_penalties[1];

        self
    }

    /// Semi-globally align a given query against the graph.
    pub fn semiglobal_banded(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, bandwidth: usize) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.x_clip,
            self.poa.y_clip,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.x_clip = MIN_SCORE;
        self.poa.y_clip = 0;

        self.query = query.to_vec();
        self.traceback = self.poa.custom_banded(query, lcsk_path, bandwidth);

        // Set the clip penalties to the original values
        self.poa.x_clip = clip_penalties[0];
        self.poa.y_clip = clip_penalties[1];
        self
    }

    /// Locally align a given query against the graph.
    pub fn local_banded(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, bandwidth: usize) -> &mut Self {
        // Store the current clip penalties
        let clip_penalties = [
            self.poa.x_clip,
            self.poa.y_clip,
        ];

        // Temporarily Over-write the clip penalties
        self.poa.x_clip = 0;
        self.poa.y_clip = 0;

        self.query = query.to_vec();
        self.traceback = self.poa.custom_banded(query, lcsk_path, bandwidth);

        // Set the clip penalties to the original values
        self.poa.x_clip = clip_penalties[0];
        self.poa.y_clip = clip_penalties[1];

        self
    }

    /// Custom align a given query against the graph with custom xclip and yclip penalties.
    pub fn custom_banded(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, bandwidth: usize) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.custom_banded(query, lcsk_path, bandwidth);
        self
    }

    pub fn custom_banded_threaded(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, bandwidth: usize, section_graph: Graph<u8, i32, Directed, usize>) -> &mut Self {
        self.poa.graph = section_graph;
        self.query = query.to_vec();
        self.traceback = self.poa.custom_banded_threaded_section(query, lcsk_path, bandwidth);
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
            let mut neighbour_nodes = self.poa.graph.neighbors_directed(node, Incoming);
            // go through the incoming neighbour nodes
            while let Some(neighbour_node) = neighbour_nodes.next() {
                let mut weight = 0;
                let neighbour_index = neighbour_node.index();
                let neighbour_score = weight_score_next_vec[neighbour_index].1;
                let mut edges = self.poa.graph.edges_connecting(neighbour_node, node);
                while let Some(edge) = edges.next() {
                    weight += edge.weight().clone();
                }
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
    pub memory_usage: usize,
}

impl Poa{

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
        Poa { match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score, x_clip: x_clip, y_clip: y_clip, graph, memory_usage: 0}
    }
    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    pub fn custom(&mut self, query: &Vec<u8>) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        // save score location of the max scoring node for the query for suffix clipping
        let mut max_in_column = vec![(0, 0); n + 1];
        let mut traceback = Traceback::with_capacity(m, n);
        let mut total_cell_usage = 0;
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
            total_cell_usage += n + 1;
            // query base and its index in the DAG (traceback matrix rows)
            for (query_index, query_base) in query.iter().enumerate() {
                let j = query_index + 1; // 0 index is initialized so we start at 1
                                         // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    let temp_score;
                    if r == *query_base {
                        temp_score = self.match_score;
                    }
                    else {
                        temp_score = self.mismatch_score;
                    }
                    TracebackCell {
                        score: traceback.get(0, j - 1).score + temp_score,
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = max(
                        TracebackCell {
                            score: MIN_SCORE,
                            op: AlignmentOperation::Match(None),
                        },
                        TracebackCell {
                            score: self.x_clip,
                            op: AlignmentOperation::Xclip(0),
                        },
                    );
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
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + temp_score,
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + self.gap_open_score,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };
                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + self.gap_open_score,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
                if max_in_column[j].0 < score.score {
                    max_in_column[j].0 = score.score;
                    max_in_column[j].1 = i;
                }
            }
        }
        // X suffix clipping
        let mut max_in_row = (0, 0);
        for j in 0..n + 1 {
            // avoid pointing to itself
            if max_in_column[j].1 == traceback.last.index() + 1 {
                continue;
            }
            let maxcell = max(
                traceback.get(traceback.last.index() + 1, j).clone(),
                TracebackCell {
                    score: max_in_column[j].0 + self.x_clip,
                    op: AlignmentOperation::Xclip(max_in_column[j].1),
                },
            );
            if max_in_row.0 < maxcell.score {
                max_in_row.0 = maxcell.score;
                max_in_row.1 = j;
            }
            traceback.set(traceback.last.index() + 1, j, maxcell);
        }
        // Y suffix clipping from the last node
        let maxcell = max(
            traceback.get(traceback.last.index() + 1, n).clone(),
            TracebackCell {
                score: max_in_row.0 + self.y_clip,
                op: AlignmentOperation::Yclip(max_in_row.1, n),
            },
        );
        if max_in_row.1 != n {
            traceback.set(traceback.last.index() + 1, n, maxcell);
        }
        self.memory_usage = (total_cell_usage * mem::size_of::<TracebackCell>()) / 1024;
        //println!("Total {}KB", (total_cell_usage * mem::size_of::<TracebackCell>()) / 1024);
        traceback
    }
    
    /// Experimental: return sequence of traversed edges
    ///
    /// Only supports alignments for sequences that have already been added,
    /// so all operations must be Match.
    pub fn edges(&self, aln: Alignment) -> Vec<usize> {
        let mut path: Vec<usize> = vec![];
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut _i: usize = 0;
        for op in aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    _i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(p);
                    let edge = self.graph.find_edge(prev, node).unwrap();
                    path.push(edge.index());
                    prev = NodeIndex::new(p);
                    _i += 1;
                }
                AlignmentOperation::Ins(None) => {}
                AlignmentOperation::Ins(Some(_)) => {}
                AlignmentOperation::Del(_) => {}
                AlignmentOperation::Xclip(_) => {}
                AlignmentOperation::Yclip(_, _) => {}
            }
        }
        path
    }

    pub fn custom_banded(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, bandwidth: usize) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        // save score location of the max scoring node for the query for suffix clipping
        let mut max_in_column = vec![(0, 0); n + 1];
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.gap_open_score, self.y_clip);
        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        let mut no_kmers = false;
        if lcsk_path.len() == 0 {
            no_kmers = true;
        }
        let mut start_banding_query_node = (0, 0);
        let mut end_banding_query_node = &(0, 0);
        let mut banding_started = false;
        let mut banding_ended = false;
        let mut current_lcsk_path_index = 0;
        let mut banded_cell_usage = 0;
        if !no_kmers {
            start_banding_query_node = lcsk_path[0];
            end_banding_query_node = lcsk_path.last().unwrap();
        }
        
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1; // 0 index is for initialization so we start at 1
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            // banding stuff do here
            let mut start = 0;
            let mut end = n;
            if !no_kmers {
                if banding_started == false {
                    //do banding till start_banding_query_node + bandwidth
                    end = start_banding_query_node.0 + bandwidth;
                }
                else if banding_ended == true {
                    // do banding till till end of table
                    start = if bandwidth > end_banding_query_node.0 {
                        0
                    } else {
                        end_banding_query_node.0 - bandwidth
                    };
                }
                else{
                    // normal banding which corrospond to path
                    // get the current lcsk path - bandwidth as start
                    start = if bandwidth > lcsk_path[current_lcsk_path_index].0 {
                        0
                    } else {
                        lcsk_path[current_lcsk_path_index].0 - bandwidth
                    };
                    // get the next lcsk path + bandwidth as end, if not available use start + 2 * bandwidth
                    if lcsk_path.len() < current_lcsk_path_index + 1 {
                        end = lcsk_path[current_lcsk_path_index + 1].0 + bandwidth;
                    }
                    else {
                        end = lcsk_path[current_lcsk_path_index].0 + bandwidth;
                    }
                }
    
                // modify the conditions banding_started, banding ended
                //println!("start {} end {} current {}", banding_started, banding_ended, current_lcsk_path_index);
                if banding_ended != true {
                    if lcsk_path[current_lcsk_path_index].1 == node.index() {
                        current_lcsk_path_index += 1;
                    }
                    if start_banding_query_node.1 == node.index() {
                        banding_started = true;
                    }
                    if end_banding_query_node.1 == node.index() {
                        banding_ended = true;
                    }
                }
            }
            //start = 0;
            if end > n {
                end = n;
            }
            //end = n;
            traceback.new_row(
                i,
                (end - start) + 1,
                self.gap_open_score,
                self.x_clip,
                start,
                end + 1,
            );
            banded_cell_usage += (end - start) + 1; 
            // query base and its index in the DAG (traceback matrix rows)
            for (query_index, query_base) in query.iter().enumerate().skip(start) {
                if query_index > end {
                    break;
                }
                let j = query_index + 1; // 0 index is initialized so we start at 1
                                         // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    let temp_score;
                    if r == *query_base {
                        temp_score = self.match_score;
                    }
                    else {
                        temp_score = self.mismatch_score;
                    }
                    TracebackCell {
                        score: traceback.get(0, j - 1).score + temp_score,
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = max(
                        TracebackCell {
                            score: MIN_SCORE,
                            op: AlignmentOperation::Match(None),
                        },
                        TracebackCell {
                            score: self.x_clip,
                            op: AlignmentOperation::Xclip(0),
                        },
                    );
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
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + temp_score,
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + self.gap_open_score,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };
                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + self.gap_open_score,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
                if max_in_column[j].0 < score.score {
                    max_in_column[j].0 = score.score;
                    max_in_column[j].1 = i;
                }
            }
        }
        // X suffix clipping
        let mut max_in_row = (0, 0);
        for j in 0..n + 1 {
            // avoid pointing to itself
            if max_in_column[j].1 == traceback.last.index() + 1 {
                continue;
            }
            let maxcell = max(
                traceback.get(traceback.last.index() + 1, j).clone(),
                TracebackCell {
                    score: max_in_column[j].0 + self.x_clip,
                    op: AlignmentOperation::Xclip(max_in_column[j].1),
                },
            );
            if max_in_row.0 < maxcell.score {
                max_in_row.0 = maxcell.score;
                max_in_row.1 = j;
            }
            traceback.set(traceback.last.index() + 1, j, maxcell);
        }
        // Y suffix clipping from the last node
        let maxcell = max(
            traceback.get(traceback.last.index() + 1, n).clone(),
            TracebackCell {
                score: max_in_row.0 + self.y_clip,
                op: AlignmentOperation::Yclip(max_in_row.1, n),
            },
        );
        if max_in_row.1 != n {
            traceback.set(traceback.last.index() + 1, n, maxcell);
        }
        // print the score here
        println!("max score is {}", traceback.get(traceback.last.index() + 1, traceback.cols).score);
        //println!("Banded {}KB", (banded_cell_usage * mem::size_of::<TracebackCell>()) / 1024);
        self.memory_usage = (banded_cell_usage * mem::size_of::<TracebackCell>()) / 1024;
        traceback
    }
    // need to convert this to use and ascending topo indices instead topo indices, extra input maybe just the hashmap
    pub fn custom_banded_threaded_section(&mut self, query: &Vec<u8>, lcsk_path: &Vec<(usize, usize)>, bandwidth: usize) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        // save score location of the max scoring node for the query for suffix clipping
        let mut max_in_column = vec![(0, 0); n + 1];
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.gap_open_score, self.y_clip);
        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        let mut no_kmers = false;
        if lcsk_path.len() == 0 {
            no_kmers = true;
        }
        let mut start_banding_query_node = (0, 0);
        let mut end_banding_query_node = &(0, 0);
        let mut banding_started = false;
        let mut banding_ended = false;
        let mut current_lcsk_path_index = 0;
        let mut banded_cell_usage = 0;
        if !no_kmers {
            start_banding_query_node = lcsk_path[0];
            end_banding_query_node = lcsk_path.last().unwrap();
        }
        //println!("Started here");
        while let Some(node) = topo.next(&self.graph) {

            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1; // 0 index is for initialization so we start at 1
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            // banding stuff do here
            let mut start = 0;
            let mut end = n;
            //println!("Start calculation");
            if !no_kmers {
                if banding_started == false {
                    //do banding till start_banding_query_node + bandwidth
                    end = start_banding_query_node.0 + bandwidth;
                }
                else if banding_ended == true {
                    // do banding till till end of table
                    start = if bandwidth > end_banding_query_node.0 {
                        0
                    } else {
                        end_banding_query_node.0 - bandwidth
                    };
                }
                else{
                    // normal banding which corrospond to path
                    // get the current lcsk path - bandwidth as start
                    start = if bandwidth > lcsk_path[current_lcsk_path_index].0 {
                        0
                    } else {
                        //println!("{}", lcsk_path[current_lcsk_path_index].0); 
                        lcsk_path[current_lcsk_path_index].0 - bandwidth
                    };
                    // get the next lcsk path + bandwidth as end, if not available use start + 2 * bandwidth
                    if lcsk_path.len() < current_lcsk_path_index + 1 {
                        //println!("{}", lcsk_path[current_lcsk_path_index + 1].0); 
                        end = lcsk_path[current_lcsk_path_index + 1].0 + bandwidth;
                    }
                    else {
                        //println!("{}", lcsk_path[current_lcsk_path_index].0); 
                        end = lcsk_path[current_lcsk_path_index].0 + bandwidth;
                    }
                }
    
                // modify the conditions banding_started, banding ended
                //println!("start {} end {} current {}", banding_started, banding_ended, current_lcsk_path_index);
                if banding_ended != true {
                    if lcsk_path[current_lcsk_path_index].1 == node.index() {
                        current_lcsk_path_index += 1;
                    }
                    if start_banding_query_node.1 == node.index() {
                        banding_started = true;
                    }
                    if end_banding_query_node.1 == node.index() {
                        banding_ended = true;
                    }
                }
            }
            //start = 0;
            if end > n {
                end = n;
            }
            //end = n;
            //println!("End calvulated!!! start {} end {} banding start {} banding end {}", start, end, banding_started, banding_ended);
            traceback.new_row(
                i,
                (end - start) + 1,
                self.gap_open_score,
                self.x_clip,
                start,
                end + 1,
            );
            banded_cell_usage += (end - start) + 1;
            //println!("End calvulated");
            // query base and its index in the DAG (traceback matrix rows)
            for (query_index, query_base) in query.iter().enumerate().skip(start) {
                if query_index > end {
                    break;
                }
                let j = query_index + 1; // 0 index is initialized so we start at 1
                                         // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    let temp_score;
                    if r == *query_base {
                        temp_score = self.match_score;
                    }
                    else {
                        temp_score = self.mismatch_score;
                    }
                    TracebackCell {
                        score: traceback.get(0, j - 1).score + temp_score,
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = max(
                        TracebackCell {
                            score: MIN_SCORE,
                            op: AlignmentOperation::Match(None),
                        },
                        TracebackCell {
                            score: self.x_clip,
                            op: AlignmentOperation::Xclip(0),
                        },
                    );
                    for prev_node in &prevs {
                        // get the prev_node index in ascending order
                        //println!("Reached here");
                        let i_p: usize = prev_node.index()+ 1; // index of previous node
                        let temp_score;
                        if r == *query_base {
                            temp_score = self.match_score;
                        }
                        else {
                            temp_score = self.mismatch_score;
                        }
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + temp_score,
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + self.gap_open_score,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                        //println!("DID not Reach here");
                    }
                    max_cell
                };
                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + self.gap_open_score,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
                //println!("DID not Reach here");
                if max_in_column[j].0 < score.score {
                    max_in_column[j].0 = score.score;
                    max_in_column[j].1 = i;
                }
            }
        }
        let mut max_in_row = (0, 0);
        for j in 0..n + 1 {
            // avoid pointing to itself
            if max_in_column[j].1 == traceback.last.index() + 1 {
                continue;
            }
            let maxcell = max(
                traceback.get(traceback.last.index() + 1, j).clone(),
                TracebackCell {
                    score: max_in_column[j].0 + self.x_clip,
                    op: AlignmentOperation::Xclip(max_in_column[j].1),
                },
            );
            if max_in_row.0 < maxcell.score {
                max_in_row.0 = maxcell.score;
                max_in_row.1 = j;
            }
            traceback.set(traceback.last.index() + 1, j, maxcell);
        }
        // Y suffix clipping from the last node
        let maxcell = max(
            traceback.get(traceback.last.index() + 1, n).clone(),
            TracebackCell {
                score: max_in_row.0 + self.y_clip,
                op: AlignmentOperation::Yclip(max_in_row.1, n),
            },
        );
        if max_in_row.1 != n {
            traceback.set(traceback.last.index() + 1, n, maxcell);
        }
        //println!("Banded {}KB", (banded_cell_usage * mem::size_of::<TracebackCell>()) / 1024);
        self.memory_usage = (banded_cell_usage * mem::size_of::<TracebackCell>()) / 1024;
        //println!("REACJED END");
        traceback
    }
    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    pub fn add_alignment(&mut self, aln: &Alignment, seq: &Vec<u8>) {
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
                        if edge_not_connected {
                            self.graph.add_edge(prev, node, 1);
                        }
                        edge_not_connected = false;
                        prev = node;
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
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
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
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                    }
                    prev = node;
                    edge_not_connected = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes and xclip
                AlignmentOperation::Xclip(_) => {}
                AlignmentOperation::Yclip(_, r) => {
                    i = *r;
                }
            }
        }
    }
}