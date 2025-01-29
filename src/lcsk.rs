use crate::bit_tree::MaxBitTree;
use fxhash::FxHasher;
use petgraph::Direction::Incoming;
use std::cmp::{min, max};
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use petgraph::Outgoing;
use petgraph::graph::NodeIndex;
use petgraph::{Directed, Graph};
use itertools::Itertools;
use fxhash::FxHashMap;
use std::simd::u8x32;

pub type POAGraph = Graph<u8, i32, Directed, usize>;
pub type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

pub fn anchoring_lcsk_path_for_threading (ascending_path: &Vec<(usize, usize)>, original_path: &Vec<(usize, usize)>, number_of_sequences: usize, graph: &POAGraph, cut_limit: usize, query_length: usize, topo_indices: Vec<usize>, query: &Vec<u8>) -> (Vec<usize>, Vec<Graph<u8, i32, Directed, usize>>, Vec<usize>, Vec<Vec<u8>>, Vec<Vec<(usize, usize)>>) {    let mut current_cut_limit = cut_limit;
    let mut section_graphs: Vec<Graph<u8, i32, Directed, usize>> = vec![];
    let mut section_queries = vec![];
    let mut section_lcsks = vec![];
    let mut temp_section_lcsk = vec![];
    let mut query_cut_off_for_lcsk;
    // add the head node to first graph
    let mut section_graph: Graph<u8, i32, Directed, usize> = Graph::default();
    let mut node_tracker: Vec<usize> = vec![0; graph.node_count()];
    let mut this_is_head_node = true;
    let mut cut_fail_times = 0;
    let cut_fail_limit = 1000; // for pacbio, to run cleanly
    let section_head_node = section_graph.add_node(graph.raw_nodes()[topo_indices[0]].weight);
    node_tracker[topo_indices[0]] = section_head_node.index();
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
            //println!("Now adding node {} node index {} lcsk node index {}", added_node.index(), topo_indices[topo_indices_index], node_index);
            let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(NodeIndex::new(topo_indices[topo_indices_index]), Incoming).collect();
            if this_is_head_node == false {
                for incoming_node in incoming_nodes {
                    let mut edges = graph.edges_connecting(incoming_node, NodeIndex::new(topo_indices[topo_indices_index]));
                    let mut incoming_weight = 0;
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                    let section_incoming_node = node_tracker[incoming_node.index()];
                    //println!("Here? {} {}", added_node.index(), section_incoming_node);
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
            //println!("CURRENT CUT LIMIT {} pos0 {} pos1 {}", current_cut_limit, pos.0, pos.1);
            // add the current node and add the edges from previous nodes
            //println!("current pos graph order {} node index {} query {}", pos.1, node_index, pos.0);
            let mut temp_num_seq = number_of_sequences;
            if cut_fail_times >= cut_fail_limit {
                temp_num_seq = number_of_sequences - 1;
            }
            if try_to_make_the_cut(graph, node_index, temp_num_seq) {
                // put the graph in vector and make new graph
                //println!("{:?}", Dot::new(&section_graph.map(|_, n| (*n) as char, |_, e| *e)));
                section_graphs.push(section_graph);
                section_graph = Graph::default();
                //println!("Cut successful");
                cut_is_made_at_this_point = true;
                this_is_head_node = true;
                // cut possible
                // if yes select as anchor
                anchors.push((pos.1, node_index, pos.0));
                section_ends.push(node_tracker[node_index]);
                // make the query stuff
                let query_start = min(anchors[anchors.len() - 1].2, anchors[anchors.len() - 2].2);
                let query_end = max(anchors[anchors.len() - 1].2, anchors[anchors.len() - 2].2);
                let section_query = query[query_start..query_end].to_vec();
                query_cut_off_for_lcsk = query_start;
                section_queries.push(section_query);
                cut_fail_times = 0;
                // increase the current cut limit by cutlimit
                current_cut_limit = pos.1 + cut_limit;
            }
            else {
                // do nothing
                //println!("Cut failed");
                cut_fail_times += 1;
                // cut not possible
                // if not go to the next node increase the current cutlimit by 1
                current_cut_limit = pos.1 + 1;
            }
        }
        // if a cut is made dont put that lcsk path to lcsk path
        if cut_is_made_at_this_point == false {
            //println!("THis is run!!!");
            // put the query point and graph point in the temp lcsk vec if the cut was not made
            temp_section_lcsk.push((original_path[index].0, node_tracker[node_index]));
        }
        else {
            //println!("ADD IS RUN");
            // if the cut is made put it in the section lcsks vec
            // subtract the cut off from all query points
            let mut filtered_section_lcsk = vec![];
            for index in 0..temp_section_lcsk.len() {
                //temp_section_lcsk[index].0 = temp_section_lcsk[index].0 - query_cut_off_for_lcsk;
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
        //println!("Now adding node {} node index {} lcsk node index {}", added_node.index(), topo_indices[topo_indices_index], node_index);
        let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(NodeIndex::new(topo_indices[topo_indices_index]), Incoming).collect();
        if this_is_head_node == false {
            for incoming_node in incoming_nodes {
                let mut edges = graph.edges_connecting(incoming_node, NodeIndex::new(topo_indices[topo_indices_index]));
                let mut incoming_weight = 0;
                while let Some(edge) = edges.next() {
                    incoming_weight += edge.weight().clone();
                }
                let section_incoming_node = node_tracker[incoming_node.index()];
                //println!("Here? {} {}", added_node.index(), section_incoming_node);
                section_graph.add_edge(NodeIndex::new(section_incoming_node), added_node, incoming_weight);
            }
        }
        else {
            this_is_head_node = false;
        }
        topo_indices_index += 1;
        
    }
    //println!("{:?}", Dot::new(&section_graph.map(|_, n| (*n) as char, |_, e| *e)));
    // final query section
    let query_start;
    if anchors.len() > 2 {
        query_start = max(anchors[anchors.len() - 1].2, anchors[anchors.len() - 2].2);
    }
    else {
        query_start = 0;
    }
    let query_end = query.len() - 1;
    let section_query = query[query_start..query_end].to_vec();
    section_queries.push(section_query);
    section_graphs.push(section_graph);
    // add the end to anchor
    anchors.push((topo_indices.len(), *topo_indices.last().unwrap(), query_length - 1));
    let mut filtered_section_lcsk = vec![];
    for index in 0..temp_section_lcsk.len() {
        //temp_section_lcsk[index].0 = temp_section_lcsk[index].0 - query_cut_off_for_lcsk;
        if temp_section_lcsk[index].0 > query_start {
            filtered_section_lcsk.push((temp_section_lcsk[index].0 - query_start, temp_section_lcsk[index].1));
        }
    }
    section_lcsks.push(filtered_section_lcsk);
    //println!("SECTION LCSK {:?}", section_lcsks);
    //println!("{:?}", Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e)));
    (section_ends, section_graphs, node_tracker, section_queries, section_lcsks)
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

pub fn lcskpp_graph(kmer_pos_vec: Vec<(u32, u32)>, kmer_path_vec: Vec<Vec<usize>>, kmers_previous_node_in_paths: Vec<Vec<u32>>, num_of_paths: usize, k: usize, kmer_graph_index: Vec<Vec<u32>>, topo_map: &Vec<usize>) -> (Vec<(usize, usize)>, Vec<(usize, usize)>, u32) {
    // return nothing if empty
    if kmer_pos_vec.is_empty() {
        return (vec![], vec![],  0);
    }

    let k = k as u32;

    let mut events: Vec<(u32, u32, u32, Vec<usize>, Vec<u32>, Vec<u32>)> = Vec::new(); // x, idx, path, prev ,kmer nodes in greph)
    let mut max_ns = vec![0; num_of_paths];
    // generate the required events
    for (idx, &(x, y)) in kmer_pos_vec.iter().enumerate() {
        // to change the code to this shit
        assert!(y == kmer_graph_index[idx][0]);
        // add the previous one as well
        // IF END SAVE y - K - 1 NODE IF START SAVE y - 1 NODE (-1 if not available)
        // added  graph node and graph node + k
        events.push((x, kmer_graph_index[idx][0], (idx + kmer_pos_vec.len()) as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone(), kmer_graph_index[idx].clone()));
        events.push((x + k - 1, *kmer_graph_index[idx].last().unwrap(), idx as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone(), kmer_graph_index[idx].clone()));
        for path in kmer_path_vec[idx].clone() {
            max_ns[path] = max(max_ns[path], x + k - 1);
            max_ns[path] = max(max_ns[path], kmer_graph_index[idx][k as usize - 1]);
        }
    }
    // ev.2 = ev.6[k as usize - 1]
    // ev.1 = ev.6[0]
    // sorting is okay with topologically converted indices
    events.sort_unstable();
    //println!("{:?}", kmer_pos_vec);
    //for event in &events {
        //println!("{:?}", event);
    //}
    //println!("DONE");
    let mut max_bit_tree_path = vec![];
    // generate empty fenwick trees
    for (_index, n) in max_ns.iter().enumerate() {
        let max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(*n as usize);
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
            //tree_update_vec = vec![];
        }
        // p is the match index
        let p = (ev.2 % kmer_pos_vec.len() as u32) as usize;
        //println!("x:{} y_start:{} y_end:{} p:{} idx:{} paths:{:?}", ev.0, ev.1, ev.2, p, ev.3, ev.4);
        // is start if higher than this
        let is_start = ev.2 >= (kmer_pos_vec.len() as u32);
        if is_start {
            //print!("IS START \n");
            dp[p].0 = k;
            dp[p].1 = -1;
            dp[p].2 = ev.5;
            dp[p].3 = ev.0;
            // go through the paths available in this event, and get the max corrosponding value and pos
            for path_index in 0..ev.3.len() {
                let path = ev.3[path_index];
                let prev_node = ev.4[path_index];
                if prev_node != u32::MAX {
                    //println!("prev node value = {}", prev_node);
                    let (temp_value, temp_position) = max_bit_tree_path[path].get(prev_node as usize);
                    //println!("temp value from fenwick tree {}", temp_value);
                    if (temp_value + k > dp[p].0) && (temp_value > 0) {
                        dp[p].0 = k + temp_value;
                        dp[p].1 = temp_position as i32;
                        best_dp = max(best_dp, (dp[p].0, p as i32, path));
                        //println!("best_dp {}", best_dp.0);
                    }
                }
            }
        } else {
            //print!("IS END \n");
            // See if this kmer continues a different kmer
            if ev.0 >= k {
                for path_index in 0..ev.3.len() {
                    let path = ev.3[path_index];
                    let prev_node = ev.4[path_index];
                    if prev_node != u32::MAX {
                        if let Ok(cont_idx) = kmer_pos_vec.binary_search(&(ev.0 - k, prev_node)) {
                            //println!("!!!!!!!!!!");
                            let prev_score = dp[cont_idx].0;
                            //let candidate = (prev_score + 1, cont_idx as i32, prev_path);
                            if prev_score + 1 > dp[p].0 {
                                //println!("CONTINUTING VALUE {}", dp[p].2.len());
                                dp[p].0 = prev_score + 1;
                                dp[p].1 = cont_idx as i32;
                                dp[p].3 = dp[p].3 + 1;
                                dp[p].2 = vec![*dp[p].2.last().unwrap()];
                            }
                            best_dp = max(best_dp, (dp[p].0, p as i32, path));
                            //println!("candidate location {} {}", ev.0 - k, prev_node);
                            //println!("cont value from candidate p {} score {}", cont_idx, dp[p].0);
                            //println!("best_dp {}", best_dp.0);
                        }
                    }
                }
            }
            // set all trees which have this match as this // maybe update this in the next iteration to prevent query overlapping
            //  update required, current x, value,  index, paths (trees)
                // ev.2 = ev.6[k as usize - 1]
    // ev.1 = ev.6[0]
            tree_update_vec.push((ev.0, (dp[p].0, p as u32), ev.5[k as usize - 1], ev.3.clone()));
            tree_update_required_level = ev.0;
        }
    }
    let mut traceback = Vec::new();
    let (best_score, mut prev_match, mut _path) = best_dp;
    //println!("BEST SCORE: {} PREV_MATCH: {}", best_score, prev_match);
    //println!("PATHS");
    
    let mut query_graph_path = vec![];
    let mut unconverted_query_graph_path = vec![];
    let mut last_node = usize::MAX;
    while prev_match >= 0 {
        //println!("{} ", prev_match);
        traceback.push(prev_match as usize);
        dp[prev_match as usize].2.reverse();
        let mut query_pos = dp[prev_match as usize].3 + dp[prev_match as usize].2.len() as u32  - 1;
        //println!("ORIGINAL Q POS {}", query_pos);
        for node in &dp[prev_match as usize].2 {
            
            let converted_node = topo_map[*node as usize];
            let current_node = *node as usize;
            if last_node == usize::MAX {
                last_node = current_node;
            }
            //println!("q pos {}", query_pos);
            query_graph_path.push((query_pos as usize, converted_node));
            unconverted_query_graph_path.push((query_pos as usize, *node as usize));
            query_pos -= 1;
            //println!("{} != {}", last_node, current_node);
            assert!(last_node >= current_node);
            last_node = current_node;
        }
        //println!("");
        prev_match = dp[prev_match as usize].1;
    }
    query_graph_path.reverse();
    unconverted_query_graph_path.reverse();
    //println!("{:?}", query_graph_path);
    (query_graph_path, unconverted_query_graph_path, best_score)
}

pub fn better_find_kmer_matches(query: &[u8], graph_sequences: &Vec<Vec<u8>>, graph_ids: &Vec<Vec<usize>>, k: usize) -> (Vec<(u32, u32)>, Vec<Vec<usize>>, Vec<Vec<u32>>, Vec<Vec<u32>>) {
    // hash the query
    let set = hash_kmers_2(query, k);
    // go through the paths and get the indices of path and make a list with query index, graph index, paths
    // aggregated result
    let mut all_result_per_path: Vec<Vec<(u32, u32, u32, u32, Vec<u32>)>> = vec![vec![]; graph_sequences.len()]; // seq, graph, graph + k, prev node in path, path index
    let mut kmers_result_vec: Vec<(u32, u32)> = vec![];
    let mut kmers_paths: Vec<Vec<usize>> = vec![];
    let mut kmers_previous_node_in_paths: Vec<Vec<u32>> = vec![];
    let mut kmer_graph_path: Vec<Vec<u32>>= vec![];
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
            let graph_index_plus_k = graph_ids[index][a_match.1 as usize - 1 + k] as u32;
            if a_match.1 > 0 {
                graph_index_minus_1 = graph_ids[index][a_match.1 as usize - 1] as u32;
            }
            else {
                graph_index_minus_1 = u32::MAX;
            }
            all_result_per_path[index].push((a_match.0, graph_index, graph_index_plus_k, graph_index_minus_1, graph_path));
        }
    }
    // using all result get the required results
    let mut loc_to_data: HashMapFx<(u32, u32), (u32, Vec<u32>, Vec<usize>, Vec<u32>)> = HashMapFx::default();
    for (index, path_result) in all_result_per_path.iter().enumerate() {
        for result_entry in path_result {
            let key = (result_entry.0, result_entry.1);
            // check if in hash map
            match loc_to_data.get_mut(&key){
                Some(x) => {
                    // if in hash map, add to prev node and path vecs
                    x.1.push(result_entry.3);
                    x.2.push(index);
                },
                None => {
                    // if not in hash map add to hash map creating prev node vec and path vec
                    loc_to_data.insert(key, (result_entry.2, vec![result_entry.3], vec![index], result_entry.4.clone()));
                },
            }
        }
    }
    for (key, value) in loc_to_data.iter().sorted().into_iter() {
        kmers_result_vec.push(*key);
        kmers_paths.push(value.2.clone());
        kmers_previous_node_in_paths.push(value.1.clone());
        kmer_graph_path.push(value.3.clone());
        //println!("{:?} / {:?}", key, value);
    }
    (kmers_result_vec, kmers_paths, kmers_previous_node_in_paths, kmer_graph_path)
}

pub fn find_sequence_in_graph (sequence: Vec<u8>, graph: &POAGraph, topo_indices: &Vec<usize>, topo_map: &HashMap<usize, usize>, error_index: usize) -> (bool, Vec<usize>, Vec<u8>) {
    let mut current_node = 0;
    // created the visit vec
    let mut visited_node: Vec<bool> = vec![false; graph.node_count() + 1];
    // created the stack to keep track
    let mut node_stack: Vec<(usize, usize)> = vec![];
    let mut current_index = 0;
    let mut final_path = vec![0; sequence.len()];
    let mut final_sequence = vec![0; sequence.len()];
    let mut error_occured = false;
    // find the first match to index 0 of sequence, if there is error in processing choose the next topo skipping error one
    let mut current_error_index = 0;
    for (_index, topo_index) in topo_indices.iter().enumerate() {
        if graph.raw_nodes()[*topo_index].weight == sequence[current_index] {
            current_node = *topo_index;
            //println!("did break at index {}", index);
            if error_index == current_error_index {
                break;
            }
            current_error_index += 1;
        }
    }
    loop {
        //println!("current {}", current_node);
        let current_node_mapped = *topo_map.get(&current_node).unwrap();
        // push to vecs required stuff
        final_path[current_index] = current_node_mapped;
        final_sequence[current_index] = graph.raw_nodes()[current_node].weight;
        // break if end of sequence reached
        if current_index >= sequence.len() - 1 {
            //println!("broke here 1");
            break;
        }
        // check if visited if not add neigbours to stack
        if !visited_node[current_node_mapped] {
            visited_node[current_node_mapped] = true;
            for neighbour in graph.neighbors(NodeIndex::new(current_node)) {
                // check if the neighbouring nodes are in sequence
                if graph.raw_nodes()[neighbour.index()].weight == sequence[current_index + 1] {
                    // add to stack the node index and current index
                    node_stack.push((neighbour.index(), current_index + 1));
                }
            }
        }
        
        // pop one from stack and process update index
        if node_stack.len() > 0  {
            (current_node, current_index) = node_stack.pop().unwrap();
        }
        // break if stack is empty
        else {
            //println!("broke here 2");
            break;
        }
    }
    if current_index != sequence.len() - 1 {
        error_occured = true;
        //println!("ERROR");
    }
    (error_occured, final_path, final_sequence)
}

pub fn dfs_get_sequence_paths (
    index: usize,
    mut temp_sequences: Vec<String>,
    graph: &POAGraph, 
    start: usize,
    current_path: Vec<usize>,
    current_sequence: Vec<u8>,
    all_paths: &mut Vec<Vec<usize>>,
    all_sequences: &mut Vec<Vec<u8>>,
    topo_map: &HashMap<usize, usize>
) {
    // go through the temp sequences and pop the ones which do not have the current base.
    let seq_len = temp_sequences.len();
    for seq_index in  0..seq_len {
        let reverse_index = seq_len - seq_index - 1;
        println!("{}", reverse_index);
        if temp_sequences[reverse_index].as_bytes()[index] != graph.raw_nodes()[start].weight {
            temp_sequences.remove(reverse_index);
            println!("removed {} ", reverse_index);
        }
    }
    if temp_sequences.len() == 0 {
        println!("removed both {}", index);
    }
    
    let mut path = current_path.clone();
    let mut sequence = current_sequence.clone();
    path.push(*topo_map.get(&start).unwrap());
    sequence.push(graph.raw_nodes()[start].weight);
    // if no neighbours, last node reached
    if (graph.neighbors(NodeIndex::new(start)).count() == 0) && (temp_sequences.len() > 0) {
        all_paths.push(path.clone());
        all_sequences.push(sequence.clone());
    } else {
        // go through neighbours recursively
        if temp_sequences.len() > 0 {
            for neighbor in graph.neighbors(NodeIndex::new(start)) {
                dfs_get_sequence_paths(index + 1,  temp_sequences.clone(), graph, neighbor.index(), path.clone(), sequence.clone(), all_paths, all_sequences, topo_map);
            }
        }
    }
}

// incomplete
pub fn lcskpp_graph_for_divided (
    kmer_pos_vec: Vec<(u32, u32, usize)>,
    kmers_plus_k: Vec<u32>,
    kmer_path_vec: Vec<Vec<(usize, usize)>>,
    kmers_previous_node_in_paths: Vec<Vec<u32>>,
    max_num_paths_per_section: usize,
    k: usize
    ) {
    // return nothing if empty
    if kmer_pos_vec.is_empty() {
        return;
    }

    let k = k as u32;

    let mut events: Vec<(u32, u32, u32, u32, Vec<(usize, usize)>, Vec<u32>)> = Vec::new(); // x, y_start, y_end, idx, (section, section_path), prev nodes
    let mut max_ns = vec![0; max_num_paths_per_section];
    // generate the required events
    for (idx, &(x, y, _section)) in kmer_pos_vec.iter().enumerate() {
        let y_plus_k = kmers_plus_k[idx];
        // add the previous one as well
        // IF END SAVE y - K - 1 NODE IF START SAVE y - 1 NODE (-1 if not available)
        events.push((x, y, y_plus_k, (idx + kmer_pos_vec.len()) as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone()));
        events.push((x + k - 1, y, y_plus_k, idx as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone()));
        for path in kmer_path_vec[idx].clone() {
            //println!("{} {}", path.1, max_num_paths_per_section);
            max_ns[path.1] = max(max_ns[path.1], x + k - 1);
            max_ns[path.1] = max(max_ns[path.1], y_plus_k);
        }
    }
    // sorting is okay with topologically converted indices
    events.sort_unstable();
    for event in events {
        println!("{:?}", event);
    }
    /* 
    let mut max_bit_tree_path = vec![];
    // generate empty fenwick trees
    for (_index, n) in max_ns.iter().enumerate() {
        let max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(*n as usize);
        max_bit_tree_path.push(max_col_dp);
    }
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (k, 0, 0); // score, coloumn, path
    // update trees only after the next query point is hit
    let mut tree_update_required_level = 0;
    let mut tree_update_vec: Vec<(u32, (u32, u32), u32, Vec<usize>)> = vec![(0, (0, 0), 0, vec![])]; //current x, value,  index, paths (trees)
    dp.resize(events.len(), (0, 0));
    for ev in events {
        if tree_update_required_level != ev.0 && tree_update_vec.len() != 0 {
            // update trees here
            for tree_info in &tree_update_vec {
                for path in &tree_info.3 {
                    max_bit_tree_path[*path].set(tree_info.2 as usize, tree_info.1);
                }
            }
            //tree_update_vec = vec![];
        }
        // p is the match index
        let p = (ev.3 % kmer_pos_vec.len() as u32) as usize;
        //println!("x:{} y_start:{} y_end:{} p:{} idx:{} paths:{:?}", ev.0, ev.1, ev.2, p, ev.3, ev.4);
        // is start if higher than this
        let is_start = ev.3 >= (kmer_pos_vec.len() as u32);
        if is_start {
            //print!("IS START \n");
            dp[p] = (k, -1);
            // go through the paths available in this event, and get the max corrosponding value and pos
            for path_index in 0..ev.4.len() {
                let path = ev.4[path_index];
                let prev_node = ev.5[path_index];
                if prev_node != u32::MAX {
                    //println!("prev node value = {}", prev_node);
                    let (temp_value, temp_position) = max_bit_tree_path[path.0].get(prev_node as usize);
                    //println!("temp value from fenwick tree {}", temp_value);
                    if (temp_value + k > dp[p].0) && (temp_value > 0) {
                        dp[p] = (k + temp_value, temp_position as i32);
                        best_dp = max(best_dp, (dp[p].0, p as i32, path.0));
                        //println!("best_dp {}", best_dp.0);
                    }
                }
            }
        } else {
            //print!("IS END \n");
            // See if this kmer continues a different kmer
            if ev.0 >= k {
                for path_index in 0..ev.4.len() {
                    let path = ev.4[path_index];
                    let prev_node = ev.5[path_index];
                    if prev_node != u32::MAX {
                        if let Ok(cont_idx) = kmer_pos_vec.binary_search(&(ev.0 - k, prev_node)) {
                            //println!("!!!!!!!!!!");
                            let prev_score = dp[cont_idx].0;
                            let candidate = (prev_score + 1, cont_idx as i32);
                            dp[p] = max(dp[p], candidate);
                            best_dp = max(best_dp, (dp[p].0, p as i32, path));
                            //println!("candidate location {} {}", ev.0 - k, prev_node);
                            //println!("cont value from candidate p {} score {}", cont_idx, candidate.0);
                            //println!("best_dp {}", best_dp.0);
                        }
                    }
                }
            }
            // set all trees which have this match as this // maybe update this in the next iteration to prevent query overlapping
            //  update required, current x, value,  index, paths (trees)
            tree_update_vec.push((ev.0, (dp[p].0, p as u32), ev.2, ev.4.clone()));
            tree_update_required_level = ev.0;
        }
    }
    let mut traceback = Vec::new();
    let (best_score, mut prev_match, mut path) = best_dp;
    //println!("BEST SCORE: {} PREV_MATCH: {}", best_score, prev_match);
    //println!("PATHS");
    while prev_match >= 0 {
        //println!("{} ", prev_match);
        traceback.push(prev_match as usize);
        prev_match = dp[prev_match as usize].1;
    }
    //traceback.reverse();
    //path :traceback;
    //score: best_score,
    //dp_vector: dp,
    best_score*/
}

pub fn find_kmer_matches_for_divided (
    query: &[u8],
    graph_sequences: &Vec<Vec<Vec<u8>>>,
    graph_ids: &Vec<Vec<Vec<usize>>>,
    k: usize
    ) -> (Vec<(u32, u32, usize)>, Vec<u32>, Vec<Vec<(usize, usize)>>, Vec<Vec<u32>>) {
    // hash the query
    let set = hash_kmers(query, k);
    println!("Hashed query");
    // go through the paths and get the indices of path and make a list with query index, graph index, paths
    let mut kmers_result_vec: Vec<(u32, u32, usize)> = vec![]; //
    let mut kmers_plus_k: Vec<u32> = vec![];
    let mut kmers_paths: Vec<Vec<(usize, usize)>> = vec![];
    let mut kmers_previous_node_in_paths: Vec<Vec<u32>> = vec![];
    for (top_index, section) in graph_sequences.iter().enumerate() {
        for (bottom_index, seq) in section.iter().enumerate() {
            let matches = find_kmer_matches_seq1_hashed(&set, seq, k);
            // go through the matches and see if they are in the already made list
            for a_match in matches {
                // first get the graph index for the match
                let graph_index = graph_ids[top_index][bottom_index][a_match.1 as usize] as u32;
                let graph_index_minus_1;
                let graph_index_plus_k = graph_ids[top_index][bottom_index][a_match.1 as usize - 1 + k] as u32;
                if a_match.1 > 0 {
                    graph_index_minus_1 = graph_ids[top_index][bottom_index][a_match.1 as usize - 1] as u32;
                }
                else {
                    graph_index_minus_1 = u32::MAX;
                }
                // if it is in the result vec, add path if indices match and path different
                if let Some(index_result) = kmers_result_vec.iter().position(|r| (r.0, r.1, r.2) == (a_match.0, graph_index, top_index)) {
                    kmers_paths[index_result].push((top_index, bottom_index));
                    kmers_previous_node_in_paths[index_result].push(graph_index_minus_1);
                }
                // not in the result vec just add
                else {
                    kmers_result_vec.push((a_match.0, graph_index, top_index));
                    kmers_plus_k.push(graph_index_plus_k);
                    kmers_paths.push(vec![(top_index, bottom_index)]);
                    kmers_previous_node_in_paths.push(vec![graph_index_minus_1]);
                }
            }
        }
    }
    
    (kmers_result_vec, kmers_plus_k, kmers_paths, kmers_previous_node_in_paths)
}

pub fn divide_poa_graph_get_paths (output_graph: &POAGraph, topo_indices: &Vec<usize>, total_num_sequences: usize, cut_threshold: usize, topo_map: &HashMap<usize, usize>) -> (Vec<Vec<Vec<usize>>>, Vec<Vec<Vec<u8>>>, usize) {
    let mut cut_start_end: Vec<(usize, usize)>= vec![];
    let mut current_topo_indices_index: usize = 0;
    let mut current_cut_limit: usize = cut_threshold;
    let mut current_start: usize = topo_indices[current_topo_indices_index];
    // add nodes to the processed nodes until we reach the cut limit
    while topo_indices.len() > current_topo_indices_index + 1 {
        if current_topo_indices_index + 1 < topo_indices.len() {
            if current_topo_indices_index >= current_cut_limit {
                println!("Try cutting...");
                let cut_result = try_to_make_the_cut(output_graph, topo_indices[current_topo_indices_index], total_num_sequences);
                // make a cut if successful update the current_cut_limit
                if cut_result {
                    println!("Cut successful! :) at {}", topo_indices[current_topo_indices_index]);
                    // this node becomes the end node make a entry in the vec
                    cut_start_end.push((current_start, topo_indices[current_topo_indices_index]));
                    current_start = topo_indices[current_topo_indices_index + 1];
                    current_cut_limit += cut_threshold;
                }
                // else just add entry to the processed nodes
                else {
                    println!("Cut failed! :( at {}", topo_indices[current_topo_indices_index]);
                    current_cut_limit += 1;
                }
            }
            println!("index {} node {} cut limit {}", current_topo_indices_index, topo_indices[current_topo_indices_index], current_cut_limit);
        }
        else {
            break;
        }
        current_topo_indices_index += 1;
    }
    cut_start_end.push((current_start, topo_indices[current_topo_indices_index]));
    println!("{:?}", cut_start_end);
    let mut all_all_paths: Vec<Vec<Vec<usize>>> = vec![];
    let mut all_all_sequences: Vec<Vec<Vec<u8>>> = vec![];
    let mut max_number_of_paths_per_section: usize = 0;
    for current_start_end in cut_start_end {
        let (start, end) = current_start_end;
        let mut all_paths: Vec<Vec<usize>> = vec![];
        let mut all_sequences: Vec<Vec<u8>> = vec![];
        println!("DOING DFS start {} end {}", start, end);
        simple_dfs_with_start_end(output_graph, start, end, vec![], vec![], &mut all_paths, &mut all_sequences, topo_map);
        println!("{:?}", all_paths);
        if max_number_of_paths_per_section < all_paths.len() {
            max_number_of_paths_per_section = all_paths.len()
        }
        all_all_paths.push(all_paths);
        all_all_sequences.push(all_sequences);
    }
    return (all_all_paths, all_all_sequences, max_number_of_paths_per_section);
}

pub fn simple_dfs_with_start_end (
    graph: &POAGraph, 
    start: usize,
    end: usize,
    current_path: Vec<usize>,
    current_sequence: Vec<u8>,
    all_paths: &mut Vec<Vec<usize>>,
    all_sequences: &mut Vec<Vec<u8>>,
    topo_map: &HashMap<usize, usize>
) {
    let mut path = current_path.clone();
    let mut sequence = current_sequence.clone();
    path.push(*topo_map.get(&start).unwrap());
    sequence.push(graph.raw_nodes()[start].weight);
    // if no neighbours, last node reached
    if (start == end) || (graph.neighbors(NodeIndex::new(start)).count() == 0) {
        // dont check for 
        if start == end {
            all_paths.push(path.clone());
            all_sequences.push(sequence.clone());
        }
    } else {
        // go through neighbours recursively
        for neighbor in graph.neighbors(NodeIndex::new(start)) {
            simple_dfs_with_start_end(graph, neighbor.index(), end, path.clone(), sequence.clone(), all_paths, all_sequences, topo_map);
        }
    }
}

pub fn find_kmer_matches(query: &[u8], graph_sequences: &Vec<Vec<u8>>, graph_ids: &Vec<Vec<usize>>, k: usize) -> (Vec<(u32, u32)>, Vec<u32>, Vec<Vec<usize>>, Vec<Vec<u32>>) {
    // hash the query
    let set = hash_kmers(query, k);
    // go through the paths and get the indices of path and make a list with query index, graph index, paths
    let mut kmers_result_vec: Vec<(u32, u32)> = vec![];
    let mut kmers_plus_k: Vec<u32> = vec![];
    let mut kmers_paths: Vec<Vec<usize>> = vec![];
    let mut kmers_previous_node_in_paths: Vec<Vec<u32>> = vec![];
    for (index, seq) in graph_sequences.iter().enumerate() {
        let matches = find_kmer_matches_seq1_hashed(&set, seq, k);
        // go through the matches and see if they are in the already made list
        for a_match in matches {
            // first get the graph index for the match
            let graph_index = graph_ids[index][a_match.1 as usize] as u32;
            let graph_index_minus_1;
            let graph_index_plus_k = graph_ids[index][a_match.1 as usize - 1 + k] as u32;
            if a_match.1 > 0 {
                graph_index_minus_1 = graph_ids[index][a_match.1 as usize - 1] as u32;
            }
            else {
                graph_index_minus_1 = u32::MAX;
            }
            // if it is in the result vec, add path if indices match and path different
            if let Some(index_result) = kmers_result_vec.iter().position(|r| (r.0, r.1) == (a_match.0, graph_index)) {
                kmers_paths[index_result].push(index);
                kmers_previous_node_in_paths[index_result].push(graph_index_minus_1);
            }
            // not in the result vec just add
            else {
                kmers_result_vec.push((a_match.0, graph_index));
                kmers_plus_k.push(graph_index_plus_k);
                kmers_paths.push(vec![index]);
                kmers_previous_node_in_paths.push(vec![graph_index_minus_1]);
            }
        }
    }
    (kmers_result_vec, kmers_plus_k, kmers_paths, kmers_previous_node_in_paths)
}

pub fn simple_dfs_all_paths (
    graph: &POAGraph, 
    start: usize,
    current_path: Vec<usize>,
    current_sequence: Vec<u8>,
    all_paths: &mut Vec<Vec<usize>>,
    all_sequences: &mut Vec<Vec<u8>>,
    topo_map: &HashMap<usize, usize>
) {
    let mut path = current_path.clone();
    let mut sequence = current_sequence.clone();
    path.push(*topo_map.get(&start).unwrap());
    sequence.push(graph.raw_nodes()[start].weight);
    // if no neighbours, last node reached
    if graph.neighbors(NodeIndex::new(start)).count() == 0 {
        all_paths.push(path.clone());
        all_sequences.push(sequence.clone());
    } else {
        // go through neighbours recursively
        for neighbor in graph.neighbors(NodeIndex::new(start)) {
            simple_dfs_all_paths(graph, neighbor.index(), path.clone(), sequence.clone(), all_paths, all_sequences, topo_map);
        }
    }
}

pub fn hash_kmers(seq: &[u8], k: usize) -> HashMapFx<&[u8], Vec<u32>> {
    let slc = seq;
    let mut set: HashMapFx<&[u8], Vec<u32>> = HashMapFx::default();
    for i in 0..(slc.len() + 1).saturating_sub(k) {
        set.entry(&slc[i..i + k])
            .or_insert_with(Vec::new)
            .push(i as u32);
    }
    set
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

pub fn find_kmer_matches_seq1_hashed(
    seq1_set: &HashMapFx<&[u8], Vec<u32>>,
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
    // no need to sort, as sorted later, test
    //matches.sort_unstable();
    matches
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
    // no need to sort, as sorted later, test
    //matches.sort_unstable();
    matches
}
