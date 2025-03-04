use crate::bit_tree::MaxBitTree;
use fxhash::FxHasher;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use petgraph::{Directed, Graph, visit::Topo, Incoming, Outgoing};
use petgraph::graph::NodeIndex;
use itertools::Itertools;
use fxhash::FxHashMap;

pub type POAGraph = Graph<u8, i32, Directed, usize>;
pub type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

pub fn threaded_lcsk_pipeline (output_graph: &POAGraph, query: &Vec<u8>, kmer_size: usize, all_paths: &Vec<Vec<usize>>, all_bases: &Vec<Vec<u8>>, cut_limit: usize) -> (Vec<(usize, usize)>, Vec<Graph<u8, i32, Directed, usize>>, Vec<Vec<usize>>, Vec<Vec<u8>>, Vec<Vec<(usize, usize)>>) {
    let mut topo = Topo::new(&output_graph);
    let mut index_topo_index_value_graph_index_vec = vec![];
    let mut key_graph_index_value_topo_index_map = HashMap::new();
    let mut incrementing_index: usize = 0;
    let mut temp_paths_converted = vec![];
    // convert the all paths to topo indices for the current iteration of graph
    while let Some(node) = topo.next(&output_graph) {
        index_topo_index_value_graph_index_vec.push(node.index());
        key_graph_index_value_topo_index_map.insert(node.index(), incrementing_index);
        incrementing_index += 1;
    }
    for path in all_paths {
        let mut temp_path = vec![];
        for node_id in path {
            temp_path.push(key_graph_index_value_topo_index_map[node_id]);
        }
        temp_paths_converted.push(temp_path);
    }
    let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = find_kmer_matches(&query, &all_bases, &temp_paths_converted, kmer_size);
    let (lcsk_path, lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, temp_paths_converted.len(), kmer_size, kmer_graph_path, &index_topo_index_value_graph_index_vec);
    let (anchors, section_graphs, section_node_trackers, section_queries, section_lcsks) = anchoring_lcsk_path_for_threading (&lcsk_path_unconverted, &lcsk_path, all_paths.len(), output_graph, cut_limit,  query.len(), index_topo_index_value_graph_index_vec, &query);
    (anchors, section_graphs, section_node_trackers, section_queries, section_lcsks)
}

pub fn lcsk_pipeline (output_graph: &POAGraph, query: &Vec<u8>, kmer_size: usize, all_paths: &Vec<Vec<usize>>, all_bases: &Vec<Vec<u8>>) -> Vec<(usize, usize)> {
    let mut topo = Topo::new(&output_graph);
    let mut index_topo_index_value_graph_index_vec = vec![];
    let mut key_graph_index_value_topo_index_map = HashMap::new();
    let mut incrementing_index: usize = 0;
    let mut temp_paths_converted = vec![];
    while let Some(node) = topo.next(&output_graph) {
        index_topo_index_value_graph_index_vec.push(node.index());
        key_graph_index_value_topo_index_map.insert(node.index(), incrementing_index);
        incrementing_index += 1;
    }
    for path in all_paths {
        let mut temp_path = vec![];
        for node_id in path {
            temp_path.push(key_graph_index_value_topo_index_map[node_id]);
        }
        temp_paths_converted.push(temp_path);
    }
    let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = find_kmer_matches(&query, &all_bases, &temp_paths_converted, kmer_size);
    let (lcsk_path, _lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, temp_paths_converted.len(), kmer_size, kmer_graph_path, &index_topo_index_value_graph_index_vec);
    lcsk_path
}

pub fn anchoring_lcsk_path_for_threading (
    ascending_path: &Vec<(usize, usize)>,
    original_path: &Vec<(usize, usize)>,
    number_of_sequences: usize,
    graph: &POAGraph,
    cut_limit: usize,
    query_length: usize,
    topo_indices: Vec<usize>,
    query: &Vec<u8>)
    -> (Vec<(usize, usize)>, Vec<Graph<u8, i32, Directed, usize>>, Vec<Vec<usize>>, Vec<Vec<u8>>, Vec<Vec<(usize, usize)>>)
    {
    let mut current_cut_limit = cut_limit;
    let mut section_graphs: Vec<Graph<u8, i32, Directed, usize>> = vec![];
    let mut section_queries = vec![];
    let mut section_lcsks = vec![];
    let mut temp_section_lcsk = vec![];
    let mut query_cut_off_for_lcsk;
    // add the head node to first graph
    let mut section_graph: Graph<u8, i32, Directed, usize> = Graph::default();
    let mut node_tracker: Vec<usize> = vec![0; graph.node_count()];
    let mut section_node_trackers: Vec<Vec<usize>> = vec![];
    let mut temp_node_tracker_for_section: Vec<usize> = vec![];
    let mut this_is_head_node = true;
    let mut cut_fail_times = 0;
    let cut_fail_limit = 1000; // for pacbio, to run cleanly
    let mut topo_indices_index = 0;
    // add head node and stuff to the anchors
    let mut anchors= vec![(0, topo_indices[0], 0)]; // order in graph, graph index, query index
    let mut section_ends = vec![];
    for (index, pos) in ascending_path.iter().enumerate() {
        let node_index = original_path[index].1;
        let mut cut_is_made_at_this_point = false;
        query_cut_off_for_lcsk = 0;
        // add nodes until we reach the pos node index to the graph
        loop  {
            //add node
            // check if this node has any incoming edges if so save them to the new graphs node
            let added_node = section_graph.add_node(graph.raw_nodes()[topo_indices[topo_indices_index]].weight);
            node_tracker[topo_indices[topo_indices_index]] = added_node.index();
            temp_node_tracker_for_section.push(topo_indices[topo_indices_index]);
            let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(NodeIndex::new(topo_indices[topo_indices_index]), Incoming).collect();
            if this_is_head_node == false {
                for incoming_node in incoming_nodes {
                    let mut edges = graph.edges_connecting(incoming_node, NodeIndex::new(topo_indices[topo_indices_index]));
                    let mut incoming_weight = 0;
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                    let section_incoming_node = node_tracker[incoming_node.index()];
                    section_graph.add_edge(NodeIndex::new(section_incoming_node), added_node, incoming_weight);
                }
            }
            else {
                this_is_head_node = false;
            }
            if topo_indices[topo_indices_index] == node_index {
                topo_indices_index += 1;
                break;
            }
            topo_indices_index += 1;
        }
        // go through the ascending path until we hit limit
        if pos.1 > current_cut_limit {
            let mut temp_num_seq = number_of_sequences;
            if cut_fail_times >= cut_fail_limit {
                temp_num_seq = number_of_sequences - 1;
            }
            if try_to_make_the_cut(graph, node_index, temp_num_seq) {
                // put the graph in vector and make new graph
                section_node_trackers.push(temp_node_tracker_for_section.clone());
                temp_node_tracker_for_section.clear();
                section_graphs.push(section_graph);
                section_graph = Graph::default();
                cut_is_made_at_this_point = true;
                this_is_head_node = true;
                // cut possible
                // if yes select as anchor
                anchors.push((pos.1, node_index, pos.0));
                // get the next node of end
                let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(NodeIndex::new(node_index), Outgoing).collect();
                assert!(incoming_nodes.len() == 1);
                section_ends.push((node_index, incoming_nodes[0].index()));
                // make the query stuff
                let mut query_start = min(anchors[anchors.len() - 1].2, anchors[anchors.len() - 2].2);
                if query_start > 0 {
                    query_start = query_start;
                }
                let query_end = max(anchors[anchors.len() - 1].2, anchors[anchors.len() - 2].2);
                let section_query = query[query_start..query_end].to_vec();
                query_cut_off_for_lcsk = query_start;
                section_queries.push(section_query);
                cut_fail_times = 0;
                // increase the current cut limit by cutlimit
                current_cut_limit = pos.1 + cut_limit;
            }
            else {
                cut_fail_times += 1;
                // cut not possible
                // if not go to the next node increase the current cutlimit by 1
                current_cut_limit = pos.1 + 1;
            }
        }
        // if a cut is made dont put that lcsk path to lcsk path
        if cut_is_made_at_this_point == false {
            // put the query point and graph point in the temp lcsk vec if the cut was not made
            temp_section_lcsk.push((original_path[index].0, node_tracker[node_index]));
        }
        else {
            // if the cut is made put it in the section lcsks vec
            // subtract the cut off from all query points
            let mut filtered_section_lcsk = vec![];
            for index in 0..temp_section_lcsk.len() {
                if temp_section_lcsk[index].0 > query_cut_off_for_lcsk {
                    filtered_section_lcsk.push((temp_section_lcsk[index].0 - query_cut_off_for_lcsk, temp_section_lcsk[index].1));
                }
            }
            section_lcsks.push(filtered_section_lcsk);
            temp_section_lcsk = vec![];
        }
    }
    // make a graph with rest of the nodes
    loop {
        if topo_indices_index == topo_indices.len() {
            break;
        }
        let added_node = section_graph.add_node(graph.raw_nodes()[topo_indices[topo_indices_index]].weight);
        node_tracker[topo_indices[topo_indices_index]] = added_node.index();
        temp_node_tracker_for_section.push(topo_indices[topo_indices_index]);
        let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(NodeIndex::new(topo_indices[topo_indices_index]), Incoming).collect();
        if this_is_head_node == false {
            for incoming_node in incoming_nodes {
                let mut edges = graph.edges_connecting(incoming_node, NodeIndex::new(topo_indices[topo_indices_index]));
                let mut incoming_weight = 0;
                while let Some(edge) = edges.next() {
                    incoming_weight += edge.weight().clone();
                }
                let section_incoming_node = node_tracker[incoming_node.index()];
                section_graph.add_edge(NodeIndex::new(section_incoming_node), added_node, incoming_weight);
            }
        }
        else {
            this_is_head_node = false;
        }
        topo_indices_index += 1;
        
    }
    // final query section
    let query_start;
    if anchors.len() > 1 {
        query_start = max(anchors[anchors.len() - 1].2, anchors[anchors.len() - 2].2);
    }
    else {
        query_start = 0;
    }
    let query_end = query.len();
    let section_query = query[query_start..query_end].to_vec();
    section_queries.push(section_query);
    section_node_trackers.push(temp_node_tracker_for_section.clone());
    section_graphs.push(section_graph);
    temp_node_tracker_for_section.clear();
    // add the end to anchor
    anchors.push((topo_indices.len(), *topo_indices.last().unwrap(), query_length - 1));
    let mut filtered_section_lcsk = vec![];
    for index in 0..temp_section_lcsk.len() {
        if temp_section_lcsk[index].0 > query_start {
            filtered_section_lcsk.push((temp_section_lcsk[index].0 - query_start, temp_section_lcsk[index].1));
        }
    }
    section_lcsks.push(filtered_section_lcsk);
    (section_ends, section_graphs, section_node_trackers, section_queries, section_lcsks)
}

pub fn try_to_make_the_cut(output_graph: &POAGraph, topo_index: usize, total_num_sequences: usize) -> bool {
    // for the selected node, see if there is only one outgoing edge
    let mut neighbour_nodes = output_graph.neighbors_directed(NodeIndex::new(topo_index), Outgoing);
    let mut number_of_outgoing_node = 0;
    let mut number_of_seq_match = false;
    while let Some(neighbour_node) = neighbour_nodes.next() {
        number_of_outgoing_node += 1;
        if number_of_outgoing_node == 1 {
            let mut edges = output_graph.edges_connecting(NodeIndex::new(topo_index), neighbour_node);
            let mut weight: i32 = 0;
            while let Some(edge) = edges.next() {
                weight += edge.weight().clone();
            }
            if weight == total_num_sequences as i32 {
                number_of_seq_match = true;
            }
        }
        if number_of_outgoing_node >= 2 {
            break;
        }
    }
    // if cut is successful return true
    if (number_of_outgoing_node == 1) && (number_of_seq_match) {
        return true;
    }
    else {
        return false;
    }
}

pub fn lcskpp_graph (
    kmer_pos_vec: Vec<(u32, u32)>,
    kmer_path_vec: Vec<Vec<usize>>,
    kmers_previous_node_in_paths: Vec<Vec<u32>>,
    num_of_paths: usize,
    k: usize,
    kmer_graph_index: Vec<Vec<u32>>,
    topo_map: &Vec<usize>)
    -> (Vec<(usize, usize)>, Vec<(usize, usize)>, u32) 
    {
    // return nothing if no kmers found
    if kmer_pos_vec.is_empty() {
        return (vec![], vec![],  0);
    }
    let k = k as u32;
    let mut events: Vec<(u32, u32, u32, Vec<usize>, Vec<u32>, Vec<u32>)> = Vec::new(); // x, idx, path, prev ,kmer nodes in greph)
    let mut max_ns = 0;
    let mut max_bit_tree_path = vec![];
    // generate the required events
    for (idx, &(x, y)) in kmer_pos_vec.iter().enumerate() {
        assert!(y == kmer_graph_index[idx][0]);
        events.push((x, kmer_graph_index[idx][0], (idx + kmer_pos_vec.len()) as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone(), kmer_graph_index[idx].clone()));
        events.push((x + k - 1, *kmer_graph_index[idx].last().unwrap(), idx as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone(), kmer_graph_index[idx].clone()));
        max_ns = max(max_ns, max(x + k - 1, *kmer_graph_index[idx].last().unwrap() + 1));
    }
    // sorting is okay with topologically converted indices
    events.sort_unstable();
    // generate empty fenwick trees, we need a tree for each path
    for _ in 0..num_of_paths {
        let max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(max_ns as usize);
        max_bit_tree_path.push(max_col_dp);
    }
    let mut dp: Vec<(u32, i32, Vec<u32>, u32)> = Vec::with_capacity(events.len()); //index is index score prev match, corrosponding graph nodes, query pos
    let mut best_dp = (k, 0, 0); // score, coloumn, path
    // update trees only after the next query point is hit
    let mut tree_update_required_level = 0;
    let mut tree_update_vec: Vec<(u32, (u32, u32), u32, Vec<usize>)> = vec![(0, (0, 0), 0, vec![])]; //current x, value,  index, paths (trees)
    dp.resize(events.len(), (0, 0, vec![], 0));
    for ev in events {
        if tree_update_required_level != ev.0 && tree_update_vec.len() != 0 {
            // update trees here
            for tree_info in &tree_update_vec {
                for path in &tree_info.3 {
                    max_bit_tree_path[*path].set(tree_info.2 as usize, tree_info.1);
                }
            }
            tree_update_vec.clear();
        }
        // p is the match index
        let p = (ev.2 % kmer_pos_vec.len() as u32) as usize;
        let is_start = ev.2 >= (kmer_pos_vec.len() as u32);
        // START
        if is_start {
            dp[p].0 = k;
            dp[p].1 = -1;
            dp[p].2 = ev.5;
            dp[p].3 = ev.0;
            // go through the paths available in this event, and get the max corrosponding value and pos
            for path_index in 0..ev.3.len() {
                let path = ev.3[path_index];
                let prev_node = ev.4[path_index];
                if prev_node != u32::MAX {
                    let (temp_value, temp_position) = max_bit_tree_path[path].get(prev_node as usize);
                    if (temp_value + k > dp[p].0) && (temp_value > 0) && (dp[p].2.first().unwrap() > dp[temp_position as usize].2.last().unwrap()) {
                        dp[p].0 = k + temp_value;
                        dp[p].1 = temp_position as i32;
                        best_dp = max(best_dp, (dp[p].0, p as i32, path));
                    }
                }
            }
        } else {
            // END
            if ev.0 >= k {
                for path_index in 0..ev.3.len() {
                    let path = ev.3[path_index];
                    let prev_node = ev.4[path_index];
                    if prev_node != u32::MAX {
                        if let Ok(cont_idx) = kmer_pos_vec.binary_search(&(ev.0 - k, prev_node)) {
                            let prev_score = dp[cont_idx].0;
                            //let candidate = (prev_score + 1, cont_idx as i32, prev_path);
                            // we should old do this if prev node is > than dp[p].2.last()
                            if (prev_score + 1 > dp[p].0) && (prev_node > *dp[p].2.last().unwrap()) {
                                dp[p].0 = prev_score + 1;
                                dp[p].1 = cont_idx as i32;
                                dp[p].3 = dp[p].3;
                                dp[p].2 = vec![*dp[p].2.last().unwrap()];
                            }
                            best_dp = max(best_dp, (dp[p].0, p as i32, path));
                        }
                    }
                }
            }
            // set all trees which have this match as this // maybe update this in the next iteration to prevent query overlapping
            //  update required, current x, value,  index, paths (trees)
            tree_update_vec.push((ev.0, (dp[p].0, p as u32), ev.5[k as usize - 1], ev.3.clone()));
            tree_update_required_level = ev.0;
        }
    }
    let (best_score, mut prev_match, mut _path) = best_dp;
    let mut query_graph_path = vec![];
    let mut unconverted_query_graph_path = vec![];
    let mut last_node = usize::MAX;
    while prev_match >= 0 {
        dp[prev_match as usize].2.reverse();
        let mut query_pos = dp[prev_match as usize].3 + k  - 1;
        for node in &dp[prev_match as usize].2 {
            let converted_node = topo_map[*node as usize];
            let current_node = *node as usize;
            if last_node == usize::MAX {
                last_node = current_node;
            }
            query_graph_path.push((query_pos as usize, converted_node));
            unconverted_query_graph_path.push((query_pos as usize, *node as usize));
            query_pos -= 1;
            assert!(last_node >= current_node);
            last_node = current_node;
        }
        prev_match = dp[prev_match as usize].1;
    }
    query_graph_path.reverse();
    unconverted_query_graph_path.reverse();
    (query_graph_path, unconverted_query_graph_path, best_score)
}

pub fn find_kmer_matches(query: &[u8], graph_sequences: &Vec<Vec<u8>>, graph_ids: &Vec<Vec<usize>>, k: usize) -> (Vec<(u32, u32)>, Vec<Vec<usize>>, Vec<Vec<u32>>, Vec<Vec<u32>>) {
    // hash the query
    let set = hash_kmers_2(query, k);
    // go through the paths and get the indices of path and make a list with query index, graph index, paths
    // aggregated result
    let mut all_result_per_path: Vec<Vec<(u32, u32, u32, Vec<u32>)>> = vec![vec![]; graph_sequences.len()]; // seq, graph, graph + k, prev node in path, path index
    let mut kmers_result_vec: Vec<(u32, u32)> = vec![];
    let mut kmers_paths: Vec<Vec<usize>> = vec![];
    let mut kmers_previous_node_in_paths: Vec<Vec<u32>> = vec![];
    let mut kmer_graph_path: Vec<Vec<u32>>= vec![];
    //println!("num of graph sequences {}", graph_sequences.len());
    for (index, seq) in graph_sequences.iter().enumerate() {
        let matches = find_kmer_matches_seq1_hashed_2(&set, seq, k);
        // go through the matches and see if they are in the already made list
        for a_match in matches {
            // get the path index and save it 
            let mut graph_path = vec![];
            // first get the graph index for the match
            for graph_node in a_match.1..a_match.1 + k as u32 {
                graph_path.push(graph_ids[index][graph_node as usize] as u32);
            }
            let graph_index = graph_ids[index][a_match.1 as usize] as u32;
            let graph_index_minus_1;
            if a_match.1 > 0 {
                graph_index_minus_1 = graph_ids[index][a_match.1 as usize - 1] as u32;
            }
            else {
                graph_index_minus_1 = u32::MAX;
            }
            all_result_per_path[index].push((a_match.0, graph_index, graph_index_minus_1, graph_path));
        }
    }
    // using all result get the required results
    let mut loc_to_data: HashMapFx<(u32, u32), (Vec<u32>, Vec<usize>, Vec<u32>)> = HashMapFx::default();
    for (index, path_result) in all_result_per_path.iter().enumerate() {
        for result_entry in path_result {
            let key = (result_entry.0, result_entry.1);
            // check if in hash map
            match loc_to_data.get_mut(&key){
                Some(x) => {
                    // if in hash map, add to prev node and path vecs
                    x.0.push(result_entry.2);
                    x.1.push(index);
                },
                None => {
                    // if not in hash map add to hash map creating prev node vec and path vec
                    loc_to_data.insert(key, (vec![result_entry.2], vec![index], result_entry.3.clone()));
                },
            }
        }
    }
    for (key, value) in loc_to_data.iter().sorted().into_iter() {
        let mut path_good = true;
        // check graph path before pushing in 
        let mut last_node = 0;
        for node in &value.2 {
            if *node > last_node {
                last_node = *node;
            }
            else {
                path_good = false;
                break;
            }
        }
        if path_good {
            kmers_result_vec.push(*key);
            kmers_paths.push(value.1.clone());
            kmers_previous_node_in_paths.push(value.0.clone());
            kmer_graph_path.push(value.2.clone());
        }
    }
    (kmers_result_vec, kmers_paths, kmers_previous_node_in_paths, kmer_graph_path)
}

pub fn hash_kmers_2(seq: &[u8], k: usize) -> FxHashMap<&[u8], Vec<u32>> {
    let mut set: FxHashMap<&[u8], Vec<u32>> = FxHashMap::default();
    
    for (i, kmer) in seq.windows(k).enumerate() {
        set.entry(kmer)
            .or_insert_with(Vec::new)
            .push(i as u32);
    }
    set
}

pub fn find_kmer_matches_seq1_hashed_2(
    seq1_set: &FxHashMap<&[u8], Vec<u32>>,
    seq2: &[u8],
    k: usize,
) -> Vec<(u32, u32)> {
    let mut matches = Vec::new();

    for i in 0..(seq2.len() + 1).saturating_sub(k) {
        let slc = &seq2[i..i + k];
        if let Some(matches1) = seq1_set.get(slc) {
            // skip non unique kmers
            if matches1.len() > 2 {
                continue;
            }
            for pos1 in matches1 {
                matches.push((*pos1, i as u32));
            }
        }
    }
    matches
}
