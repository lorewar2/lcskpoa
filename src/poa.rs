
use std::cmp::{max, Ordering};
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
        if !(self.matrix[i].1 > j || self.matrix[i].2 <= j || self.matrix[i].0.is_empty()) {
            let real_position = j - self.matrix[i].1;
            &self.matrix[i].0[real_position]
        }
        // behind the band, met the edge
        else if j == 0 {
            &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Del(None),
            }
        }
        // infront of the band
        else if j >= self.matrix[i].2 {
            &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Ins(None),
            }
        }
        // behind the band
        else {
            &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            }
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
            ops.push(self.get(i, j).op);
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
    poa: Poa,
}

impl Aligner {
    /// Create new instance.
    pub fn new(match_score: i32, mismatch_score: i32, gap_open_score: i32, reference: &Vec<u8>, x_clip: i32, y_clip: i32, _band_size: i32) -> Self {
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
    pub memory_usage: usize,
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
        //println!("Total {}KB", (total_cell_usage * mem::size_of::<TracebackCell>()) / 1024);
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