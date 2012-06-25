"""
    Copyright (C) 2009  Emre Guney, Javier Garcia-Garcia

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#########################################################################
# Graph Utility Methods
# Methods to
# - create a network
# - filter a network based on degree
# - find paths and (connected) components of given network
# - randomize a network
# - read a network from SIF file
# - analyze network 
#
#########################################################################

import networkx
import random
import copy

#MIN_NUMBER_OF_PERTURBATION = 25
MAX_NUMBER_OF_TRIAL = 6

def main():
    g = create_network_from_sif_file("interactions.sif")
    degrees = g.degree(with_labels=True)
    node_to_values = get_node_degree_related_values(g, ["v2","v3"])
    for v in g.nodes():
	print v, degrees[v], node_to_values[v]
    create_R_analyze_network_script(g, seeds = ["v2","v3"], out_path="./")
    g2=prepare_data.network_utilities.prune_graph_at_given_percentage(g,10)
    return

def create_graph_with_same_type(G):
    return create_empty_copy(G)


def create_graph():
    """
        Creates & returns a graph
    """
    return networkx.Graph()


def get_number_of_distinct_edges(G):
    edge_list = G.edges()
    edge_set = set()
    for id1, id2, data in edge_list:
        edge_set.add((id1, id2))
    return len(edge_set)


def get_shortest_path_between(G, source_id, target_id):
    return networkx.shortest_path(G, source_id, target_id)


def get_all_paths_from(G, source_id): 
    """
        get all paths from source node to all possible nodes in a dictionary
    """
    return networkx.single_source_dijkstra_path(G, source_id)


def get_path_network(G, listNodes, path_length_cutoff=10000):
    """
    Returns a subgraph containing only nodes in listNodes

    Nodes are connected if exist a path between them, and weight of the edge consists in the length of shortest path

    If the shortest between two nodes includes another node in the list, this edge is not added
    """
    # First, check all nodes in listNodes are in network
    new_graph = networkx.MultiGraph()

    for x in xrange(len(listNodes)):
        for y in xrange(x):
            sp = networkx.shortest_path(G, listNodes[x], listNodes[y])
            if sp:
                if (len(sp)-1)<=path_length_cutoff:
                    if len(set(listNodes).intersection(set(sp[1:-1])))==0:
                        new_graph.add_edge(listNodes[x],listNodes[y],len(sp)-1)
            
    return new_graph


def get_subgraph(G, nodes):
    """
	NetworkX subgraph method wrapper
    """
    return G.subgraph(nodes)


def get_connected_components(G, return_as_graph_list=True):
    """
        Finds (strongly in the case of directed network) connected components of graph
        returnAsGraphList: returns list of graph objects corresponding to connected components (from larger to smaller)
        otherwise returns list of node list corresponding nodes in connected components
    """
    result_list = []

    if return_as_graph_list:
        result_list = networkx.connected_component_subgraphs(G)
    else:
        result_list = networkx.connected_components(G)

    return result_list


def create_network_from_sif_file(network_file_in_sif, use_edge_data = False, delim = None, include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type = use_edge_data, delim = delim)
    g=networkx.Graph()
    if include_unconnected:
	g.add_nodes_from(setNode)
    if use_edge_data:
	for e,w in dictEdge.iteritems():
	    u,v = e
	    g.add_edge(u,v,w)
    else:
	g.add_edges_from(setEdge)
    return g


def output_network_in_sif(g, output_file_name, delim = " ", include_unconnected=True):
    f = open(output_file_name, 'w')
    included_nodes = set()
    for u,v in g.edges_iter():
	f.write("%s%s%s%s%s\n" % (u, delim, g.get_edge(u,v), delim, v) )
	included_nodes.add(u)
	included_nodes.add(v)
    if include_unconnected:
	for u in g.nodes_iter():
	    if u not in included_nodes:
		f.write("%s\n" % u )
    f.close()
    return


def analyze_network(g, out_file = None, seeds = None, calculate_radius = False):
    if out_file is None:
	from sys import stdout 
	out_method = stdout.write
    else:
	file = open(out_file, "a")
	out_method = file.write
    # print networkx.info(g) # currently buggy but fixed in the next version of networkx
    out_method("(V,E): %d %d\n" % (g.number_of_nodes(), g.number_of_edges()))
    out_method("V/E: %f\n" % (int(g.number_of_nodes()) / float(g.number_of_edges())))
    degrees = g.degree().values()
    degrees.sort()
    out_method("Average degree: %f\n" % float(sum(degrees)/float(len(degrees))))
    out_method("Most connected 20 nodes: %s\n" % degrees[-20:])
    connected_components = networkx.connected_components(g)
    out_method("Connected component sizes: %s\n" % map(len, connected_components))
    # Radius calculation is very time consuming
    if calculate_radius:
	out_method("Radius (of the largest connected component): %d\n" % get_network_radius(g.subgraph(connected_components[0])))
    #print get_network_degree_histogram(g)
    if seeds is not None:
	if len(seeds) <= 1:
	    out_method("Average linker degree: %f\n" % (0))
	    out_method("Average seed connecting shortest paths length: %f\n" % (0))
	else:
	    node_to_ld = get_node_linker_degrees(g, seeds)
	    vals = [ node_to_ld[v] for v in seeds ]
	    out_method("Average linker degree: %f\n" % (float(sum(vals))/len(vals)))
	    #sp_lengths = networkx.all_pairs_shortest_path_length(g)
	    sum_length = 0
	    count = 0
	    for i, u in enumerate(seeds):
		sp_lengths = networkx.single_source_dijkstra_path_length(g, u)
		for j, v in enumerate(seeds):
		    if j<i:
			sum_length += sp_lengths[v] #[u][v]
			count += 1
	    out_method("Average seed connecting shortest paths length: %f\n" % (float(sum_length)/count))
    if out_file is not None:
	file.close()
    return


def get_node_linker_degrees(g, seeds):
    """
    Get node linker degrees
    """
    node_to_ld = {}
    for current_node in g.nodes():
	value = 0 
	for current_neighbor in g.neighbors(current_node):
	    if current_neighbor!=current_node:
		if current_neighbor in seeds:
		    value += 1 
	node_to_ld[current_node] = value
    return node_to_ld


def get_node_degree_related_values(g, seeds):
    """
    Get node degree, linker degree and linker degrees at 2 level neighborhood
    Returns a dictionary of nodes to (d, ld, d2, ld2)
    """
    node_to_values = {}
    seeds = set(seeds)
    for current_node in g.nodes():
	neighbors = set(g.neighbors(current_node))
	neighbors.discard(current_node)
	neighbors2 = copy.deepcopy(neighbors)
	d = len(neighbors)
	ld = len(neighbors & seeds)
	for current_neighbor in neighbors:
	    neighbors2 |= set(g.neighbors(current_neighbor))
	if len(neighbors2) > 0: # if node is not unconnected
	    neighbors2.remove(current_node)
	d2 = len(neighbors2)
	ld2 = len(neighbors2 & seeds)
	node_to_values[current_node] = (d, ld, d2, ld2)
    return node_to_values


def filter_network(g, degree_threshold=None, largest_connected_component=True):
    print "V,E:", g.number_of_nodes(), g.number_of_edges()
    degrees = g.degree(with_labels=True)
    subgraph_nodes = []
    for id, d in degrees.iteritems():
	if degree_threshold is None or d <= degree_threshold:
	    subgraph_nodes.append(id)
    g_filtered = g.subgraph(subgraph_nodes)
    if largest_connected_component:
	component_nodes = networkx.connected_components(g_filtered)[0]
	g_filtered = g_filtered.subgraph(component_nodes )
    print "V,E filtered:", g_filtered.number_of_nodes(), g_filtered.number_of_edges()
    return g_filtered


def get_edge_values_from_sif_attribute_file(file_name, store_edge_type=False, delim=None):
    """
	store_edge_type: if True returns dict in [(u, "pp", v)] = val format, if False returns dict in [(u,v)] = val format
    """
    edge_to_values = {}
    f=open(file_name)
    line = f.readline() # read attribute name
    line = f.readline()
    while line:
	if delim is None:
	    words = line[:-1].split()
	else:
	    words = line[:-1].split(delim)
	if len(words) != 5 or words[3] != "=":
	    print "format error", line
	    continue
	if store_edge_type:
	    edge_to_values.setdefault((words[0], words[1], words[2]), set()).add(words[4])
	else:
	    edge_to_values.setdefault((words[0], words[2]), set()).add(words[4])
	line = f.readline()
    f.close()
    return edge_to_values


def get_nodes_and_edges_from_sif_file(file_name, store_edge_type = False, delim=None, data_to_float=True):
    """
	Parse sif file into node and edge sets and dictionaries
	returns setNode, setEdge, dictNode, dictEdge
	store_edge_type: if True, dictEdge[(u,v)] = edge_value
	delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    f=open(file_name)
    for line in f:
	if delim is None:
	    words = line[:-1].split()
	else:
	    words = line[:-1].split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
	    if data_to_float:
		score = float(words[1])
	    else:
		score = words[1]
            dictNode[id1] = score
        elif len(words) == 3: 
            id2 = words[2]
            setNode.add(id2)
	    setEdge.add((id1, id2))
            if store_edge_type:
		if data_to_float:
		    dictEdge[(id1, id2)] = float(words[1])
		else:
		    dictEdge[(id1, id2)] = words[1]
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    return setNode, setEdge, dictNode, dictEdge


def get_jaccard_index_map(g):
    edge_to_jaccard = {}
    for u,v in g.edges_iter():
	u_neighbors = set(g.neighbors(u))
	v_neighbors = set(g.neighbors(v))
	edge_to_jaccard[(u,v)] = float(len(u_neighbors & v_neighbors)) / len(u_neighbors | v_neighbors)
    return edge_to_jaccard 


def get_clustering_coefficient_map(g):
    return networkx.clustering(g, with_labels=True)


def get_network_radius(g):
    return networkx.radius(g)


def get_network_degree_histogram(g):
    return networkx.degree_histogram(g)


def randomize_graph(graph, randomization_type, allow_self_edges = True):
    """
    Creates a random network from given network as a networkx graph
    randomization_type: 
        - "random": add same number of edges randomly between nodes of original graph
        - "preserve_topology": keep edges, shuffle nodes of original graph
        - "preserve_topology_and_node_degree": keep edges, shuffle nodes of original graph with the nodes of same degree
        - "preserve_degree_distribution": remove an edge between two random nodes with degrees k, l then add to two nodes with degrees k-1 & l-1, then shuffle nodes
        - "preserve_degree_distribution_and_node_degree": remove 2 random edges between a-b and c-d where degree(a)=degree(c) and degree(b)=degree(d) then add 2 edges between a-d and b-c, then shuffle nodes with the same degree
	- "erdos_renyi": creates a graph where edges are redistributed based on erdos renyi random model. 
	- "barabasi_albert": creates a graph where edges are redistributed based on barabasi albert model (preferential attachment). 
    """

    debug = False

    n_node = graph.number_of_nodes()
    n_edge = graph.number_of_edges()

    if randomization_type == "erdos_renyi":
	#raise Exception("Work in progress")
	p = float(2 * n_edge) / (n_node*n_node - 2*n_node)
	# Chooses each of the possible [n(n-1)]/2 edges with probability p 
	new_graph = networkx.erdos_renyi_graph(n_node, p)
	mapping = dict(zip(new_graph.nodes(), graph.nodes()))
	new_graph = networkx.relabel_nodes(new_graph, mapping)
	available_edges = graph.edges()
	
	# Map graph from random model to new graph
        for edge in new_graph.edges():
	    if len(available_edges) > 0:
		edge_org = available_edges.pop()
		if debug:
		    print "From random:", (edge[0], edge[1])
		new_graph.add_edge(edge[0], edge[1], graph.get_edge(edge_org[0], edge_org[1]))
	    # If the random model added too many edges
	    else:
		if debug:
		    print "Removing:", edge
		new_graph.remove_edge(edge[0], edge[1])

	# If the random model failed to add enough edges
	nodes = new_graph.nodes()
	for edge_org in available_edges:
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
	    if debug:
		print "Adding:", (source_id, target_id)
	    new_graph.add_edge(source_id, target_id, graph.get_edge(edge_org[0], edge_org[1]))
	return new_graph

    if randomization_type == "barabasi_albert":
	#raise Exception("Work in progress")
	if n_edge >= n_node:
	    # A graph of n nodes is grown by attaching new nodes each with m edges that are preferentially attached to existing nodes with high degree
	    new_graph = networkx.barabasi_albert_graph(n_node, n_edge / n_node)
	    mapping = dict(zip(new_graph.nodes(), graph.nodes()))
	    new_graph = networkx.relabel_nodes(new_graph, mapping)
	else:
	    new_graph = networkx.create_empty_copy(graph) 

	available_edges = graph.edges() 
	degree_map = networkx.degree(new_graph, with_labels=True)
	nodes = new_graph.nodes()

	# Map graph from random model to new graph
        for edge in new_graph.edges():
	    if len(available_edges) > 0:
		edge_org = available_edges.pop()
		if debug:
		    print "From random:", (edge[0], edge[1])
		new_graph.add_edge(edge[0], edge[1], graph.get_edge(edge_org[0], edge_org[1]))
	    # If the random model added too many edges
	    else:
		nodes_to_select = [ id for id, d in degree_map.items() for j in xrange(d+1) ]
		source_id = random.choice(nodes())
		target_id = random.choice(nodes_to_select)
		if debug:
		    print "Removing:", (source_id, target_id)
		new_graph.remove_edge(source_id, target_id)
		degree_map[source_id] -= 1 
		degree_map[target_id] -= 1 

	# If the random model failed to add enough edges
	for edge_org in available_edges:
	    nodes_to_select = [ id for id, d in degree_map.items() for j in xrange(d+1) ]
            source_id = random.choice(nodes)
            target_id = random.choice(nodes_to_select)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
		source_id = random.choice(nodes)
		target_id = random.choice(nodes_to_select)
	    if debug:
		print "Adding:", (source_id, target_id)
	    new_graph.add_edge(source_id, target_id, graph.get_edge(edge_org[0], edge_org[1]))
	    degree_map[source_id] += 1 
	    degree_map[target_id] += 1 

	return new_graph

    new_graph = networkx.create_empty_copy(graph) 
    #new_graph.add_nodes_from(graph.nodes())

    if randomization_type == "random":
	nodes = new_graph.nodes()
        for edge in graph.edges():
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            new_graph.add_edge(source_id, target_id, graph.get_edge(edge[0],edge[1]))
        
    elif randomization_type=="preserve_topology": # shuffle_nodes
        nodes = graph.nodes()
        random_nodes = graph.nodes()
        random.shuffle(random_nodes)
        equivalences = dict([(nodes[i],random_nodes[i]) for i in xrange(len(nodes))])
        new_graph.add_edges_from([ (equivalences[current_edge[0]],equivalences[current_edge[1]],graph.get_edge(current_edge[0],current_edge[1])) for current_edge in graph.edges() ])

    elif randomization_type=="preserve_topology_and_node_degree": # shuffle_nodes_within_same_degree
        nodes_by_degree = dict( (degree,[]) for degree in graph.degree() )
        graph_degree = graph.degree(with_labels=True)
        [ nodes_by_degree[graph_degree[node]].append(node) for node in graph_degree ]
        equivalences = {}
        for current_degree in nodes_by_degree.keys():
            nodes = nodes_by_degree[current_degree]
            random_nodes = list(nodes)
            random.shuffle(random_nodes)
            equivalences.update(dict([(nodes[i],random_nodes[i]) for i in xrange(len(nodes))]))
        new_graph.add_edges_from([ (equivalences[current_edge[0]],equivalences[current_edge[1]], graph.get_edge(current_edge[0],current_edge[1])) for current_edge in graph.edges() ])
        
    elif randomization_type=="preserve_degree_distribution":
        ## add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge(current_node1, current_node2))
        max_degree = sorted(graph.degree())[-1]
        #nodes_by_degree = dict( (degree,{}) for degree in graph.degree() )
        nodes_by_degree = dict( (degree,{}) for degree in xrange(max_degree+1) )
        graph_degree = graph.degree(with_labels=True)
        [ nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree ]
        #print new_graph.nodes(), new_graph.edges()
        #print nodes_by_degree
        #if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(n_edge/2, n_edge)
        #else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)
        n_perturbation = random.randint(n_edge/2, n_edge)
        for i in xrange(n_perturbation):
            n_trial = 0
            while True:
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
		    if debug:
			print "Warning: Max number of trials exceeded in perturbation ", i
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                while source_degree < 1:  
                    source_id = random.choice(new_graph.nodes())
                    source_degree = new_graph.degree(source_id)
                target_id = random.choice(new_graph.neighbors(source_id))
                target_degree = new_graph.degree(target_id)
                del nodes_by_degree[source_degree][source_id] 
                nodes_by_degree[source_degree-1].setdefault(source_id)
                del nodes_by_degree[target_degree][target_id] 
                nodes_by_degree[target_degree-1].setdefault(target_id)
                ## not very important to check for cases where new_source = source (v.v. for targets) 
                new_source_id = random.choice(nodes_by_degree[source_degree-1].keys())
                new_target_id = random.choice(nodes_by_degree[target_degree-1].keys())
		if debug:
		    print source_id, target_id, " / ", new_source_id, new_target_id
                ## check if going to add an existing edge or self edge
                if new_graph.has_edge(new_source_id, new_target_id) or new_source_id == new_target_id:
                    del nodes_by_degree[source_degree-1][source_id] 
                    nodes_by_degree[source_degree].setdefault(source_id)
                    del nodes_by_degree[target_degree-1][target_id] 
                    nodes_by_degree[target_degree].setdefault(target_id)
                    continue
		if debug:
		    print "rm %d %d" % (source_id, target_id)
                edge_data = new_graph.get_edge(source_id, target_id)
                new_graph.delete_edge(source_id, target_id)
		if debug:
		    print "add %d %d" % (new_source_id, new_target_id)
                new_graph.add_edge(new_source_id, new_target_id, edge_data)
                del nodes_by_degree[source_degree-1][new_source_id] 
                nodes_by_degree[source_degree].setdefault(new_source_id)
                del nodes_by_degree[target_degree-1][new_target_id] 
                nodes_by_degree[target_degree].setdefault(new_target_id)
                break
        #self.randomize_graph(new_graph, "preserve_topology")

    elif randomization_type=="preserve_degree_distribution_and_node_degree":
        ## add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge(current_node1, current_node2))
        nodes_by_degree = dict( (degree,{}) for degree in graph.degree() )
        graph_degree = graph.degree(with_labels=True)
        [ nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree ]
        
        #if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(1, n_edge)
        #else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)
        n_perturbation = random.randint(n_edge/2, n_edge)
        for i in xrange(n_perturbation):
            source_id = random.choice(new_graph.nodes())
            source_degree = new_graph.degree(source_id)
            ## find a node for which another node with the same degree exists
            #available_neighbors = []
            n_trial = 0
            while True: #(len(nodes_by_degree[source_degree]) < 2 or len(available_neighbors) < 1):
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
		    if debug:
			print "Warning: Max number of trials exceeded in perturbation ", i
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                if len(nodes_by_degree[source_degree]) < 2:
                    continue
                available_neighbors = []
                ## find a neighbor for which another node with the same degree exists
                for neighbor_id in new_graph.neighbors_iter(source_id):
                    if source_degree == new_graph.degree(neighbor_id): 
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 2:
                            available_neighbors.append(neighbor_id)
                    else:
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 1:
                            available_neighbors.append(neighbor_id)
                if len(available_neighbors) < 1:
                    continue
                target_id = random.choice(available_neighbors)
                target_degree = new_graph.degree(target_id)
                ## select a new source node with different id
		n_trial2 = 0
		inner_break = False
                while True:
		    n_trial2 += 1
		    if n_trial2 > MAX_NUMBER_OF_TRIAL:
			if debug:
			    print "Warning: Max number of trials exceeded in perturbation ", i
			inner_break = True
			break
                    new_source_id = random.choice(nodes_by_degree[source_degree].keys())
                    while new_source_id == source_id:
                        new_source_id = random.choice(nodes_by_degree[source_degree].keys())
                    new_available_neighbors = []
                    ## find a neighbor as new target node for which id is different from target and has an id equivalent to target
                    for neighbor_id in new_graph.neighbors_iter(new_source_id):
                        if target_degree == new_graph.degree(neighbor_id): 
                            new_available_neighbors.append(neighbor_id)
                    if len(new_available_neighbors) < 1:
                        continue
                    new_target_id = random.choice(new_available_neighbors)
                    if len(new_available_neighbors) > 1:
                        while new_target_id == target_id:
                            new_target_id = random.choice(new_available_neighbors)
                            #print new_available_neighbors, new_target_id
                    else:
                        new_target_id = new_available_neighbors[0]
                    break
		if inner_break:
		    break
		if debug:
		    print source_id, target_id, " / ", new_source_id, new_target_id
                if source_id == new_target_id or new_source_id == target_id:
                    continue
                if new_graph.has_edge(source_id, new_target_id) or new_graph.has_edge(new_source_id, target_id):
                    continue
		if debug:
		    print "rm %d %d" % (source_id, target_id)
		    print "rm %d %d" % (new_source_id, new_target_id)
                edge_data_1 = new_graph.get_edge(source_id, target_id)
                edge_data_2 = new_graph.get_edge(new_source_id, new_target_id)
                new_graph.delete_edge(source_id, target_id)
                new_graph.delete_edge(new_source_id, new_target_id)
		if debug:
		    print "add %d %d" % (source_id, new_target_id)
		    print "add %d %d" % (new_source_id, target_id)
                new_graph.add_edge(source_id, new_target_id, edge_data_1)
                new_graph.add_edge(new_source_id, target_id, edge_data_2)

    else:
        raise Exception("Unknown randomization type %s" % randomization_type)

    return new_graph


def permute_graph_at_given_percentage(graph, percentage, allow_self_edges = True):
    """
    Randomly selects percentage% of edges and reconnects them
    """
    new_graph = graph.copy() 
    nodes = new_graph.nodes()
    edges = graph.edges()
    count = int(round(len(edges) * float(percentage) / 100))
    random.shuffle(edges)
    for edge in edges[:count]:
	new_graph.remove_edge(edge[0], edge[1])
	source_id = random.choice(nodes)
	target_id = random.choice(nodes)
	while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id) \
		    or (edge[0] == source_id and edge[1] == target_id) or (edge[1] == source_id and edge[0] == target_id):
	    source_id = random.choice(nodes)
	    target_id = random.choice(nodes)
	new_graph.add_edge(source_id, target_id, graph.get_edge(edge[0],edge[1]))
    return new_graph 
 

def prune_non_seed_interactions_at_given_percentage(graph, percentage, reserved_nodes):
    """
    Randomly selects percentage% of edges and removes them (provided that they dont
    have any connection with a node in reserved_nodes)
    """
    new_graph = graph.copy()
    nodes = new_graph.nodes()
    candidate_nodes = set(nodes) - set(reserved_nodes)
    edges = graph.edges()
    candidate_edges = [ edge for edge in edges if edge[0] in candidate_nodes and edge[1] in candidate_nodes ]
    count = int(round(len(edges) * float(percentage) / 100))
    random.shuffle(candidate_edges)
    #i = 0
    #for edge in candidate_edges:
    for edge in candidate_edges[:count]:
	#if i < count:
	#    new_graph.remove_edge(edge[0], edge[1])
	#i += 1
	new_graph.remove_edge(edge[0], edge[1])
    #if i < count:
    if len(candidate_edges) < count:
	print "Warning: Pruning percentage is not achieved due to reserved nodes"
    return new_graph 


def prune_graph_at_given_percentage(graph, percentage):
    """
    Randomly selects percentage% of edges and removes them 
    """
    new_graph = graph.copy()
    edges = graph.edges()
    count = int(round(len(edges) * float(percentage) / 100))
    random.shuffle(edges)
    for edge in edges[:count]:
	new_graph.remove_edge(edge[0], edge[1])
    return new_graph 


def create_R_analyze_network_script(g, seeds = None, out_path = "./", title = "", scale_by_log=False):
    if scale_by_log:
	f = open(out_path + "analyze_network_log_scaled.r", "w")
	f.write("postscript(\"%sanalyze_network_log_scaled.eps\", width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = \"special\", title = \"%s\")\n" % (out_path, title))
    else:
	f = open(out_path + "analyze_network.r", "w")
	f.write("postscript(\"%sanalyze_network.eps\", width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = \"special\", title = \"%s\")\n" % (out_path, title))
    f.write("par(mfrow=c(3,2))\n")
    
    node_to_values = get_node_degree_related_values(g, seeds)

    f.write("f<-function(l) { \nr<-c()\nfor(i in 1:length(l)) {\n if(l[i] <= 1) {\n  r[i]<-l[i]\n }\n else { \n  r[i]<-(1+log10(l[i]))\n } \n}\nreturn(r)\n}\n")
    # Degree distribution
    #degrees = g.degree(with_labels = True) 
    #n_node = len(degrees)
    #d_max = max(degrees.values())
    d_max = max(node_to_values.values(), key=lambda x: x[0])[0]
    non_seed_degree_counts = [ 0 ] * (d_max + 1)
    seed_degree_counts = [ 0 ] * (d_max + 1)
    #for v,d in degrees.iteritems():
    for v, values in node_to_values.iteritems():
	d, ld, d2, ld2 = values
	#if d_ != d: print v, d, d_ # networkx includes self edges in degree
	if v in seeds:
	    seed_degree_counts[d] += 1
	else:
	    non_seed_degree_counts[d] += 1
    f.write("n<-c(%s)\n" % ",".join(map(str, seed_degree_counts)))
    f.write("N<-c(%s)\n" % ",".join(map(str, non_seed_degree_counts)))
    f.write("s<-(n+N)\n")
    if scale_by_log:
	f.write("A<-matrix(c(f(n), f(N+n)-f(n)), nrow=2, byrow=TRUE)\n")
	f.write("barplot(A, yaxt=\"n\", names.arg=c(0:(length(n)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Degree\", ylab=\"Number of nodes (Log scale)\")\n")
	f.write("axis(2, 0:5, labels=c(0,%s), las=2)\n" % ",".join(map(str, [ 10**i for i in xrange(5) ])))
    else:
	f.write("A<-matrix(c(n, N), nrow=2, byrow=TRUE)\n")
	##f.write("barplot(A, yaxp=c(0,max(N+n),1), xaxp=c(1, length(n), 1), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Degree\", ylab=\"Number of nodes\")\n")
	f.write("barplot(A, ylim=c(0,round(max(s)/4)), names.arg=c(0:(length(n)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Degree\", ylab=\"Number of nodes\")\n")
	#f.write("barplot(A, xaxt=\"n\", legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Degree\", ylab=\"Number of nodes\")\n")
	#f.write("axis(1, seq(0,length(n)-1,round(length(n)/5)), labels=seq(0,length(n)-1,round(length(n)/5)))\n")

    # Log-log plot of degree distribution
    #f.write("par(fig=c(0.6, 1, 0.1, 0.5), new = T)\n") # Make this plot inner
    f.write("p<-c()\n")
    f.write("d<-c()\n")
    f.write("for(i in 2:length(s)) { if(s[i] != 0) { p[i]<-s[i]; d[i]<-i; } }\n")
    f.write("fit<-lm(log10(d) ~ log10(p))\n")
    f.write("plot(0:(length(s)-1), s, log=\"xy\", yaxt=\"n\", xlab=\"Degree (Log scale)\", ylab=\"Number of nodes (Log scale)\")\n")
    f.write("axis(2, las=2)\n")
    f.write("#plot(d, p, yaxt=\"n\", xlab=\"log(degree)\", ylab=\"log(# of nodes)\")\n")
    f.write("#axis(2, 0:5, labels=c(0,%s), las=2)\n" % ",".join(map(str, [ 10**i for i in xrange(5) ])))
    f.write("abline(fit, col=2, lty=2)\n")
    f.write("legend(\"topright\", c(\"linear (ls) fit\"), col=c(2), lty=2)\n")
    #f.write("par(fig=c(0, 1, 0, 1))\n")

    # Linker degree distribution
    #node_to_ld = get_node_linker_degrees(g, seeds)
    #ld_max = max(node_to_ld.values())
    ld_max = max(node_to_values.values(), key=lambda x: x[1])[1]
    non_seed_ldegree_counts = [ 0 ] * (ld_max + 1)
    seed_ldegree_counts = [ 0 ] * (ld_max + 1)
    #for v,ld in node_to_ld.iteritems():
    for v, values in node_to_values.iteritems():
	d, ld, d2, ld2 = values
	#d, ld_, d2, ld2 = node_to_values[v]
	#assert (ld_ == ld)
	#if ld_ != ld: print v, ld, ld_
	if v in seeds:
	    seed_ldegree_counts[ld] += 1
	else:
	    non_seed_ldegree_counts[ld] += 1
    f.write("ln<-c(%s)\n" % ",".join(map(str, seed_ldegree_counts)))
    f.write("lN<-c(%s)\n" % ",".join(map(str, non_seed_ldegree_counts)))
    if scale_by_log:
	f.write("lA<-matrix(c(f(ln), f(lN+ln)-f(ln)), nrow=2, byrow=TRUE)\n")
	f.write("barplot(lA, yaxt=\"n\", names.arg=c(0:(length(ln)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Linker Degree\", ylab=\"Number of nodes (Log scale)\")\n")
	f.write("#lA<-matrix(c(f(lN+ln), f(ln)), nrow=2, byrow=TRUE)\n")
	f.write("#barplot(lA, beside=TRUE, yaxt=\"n\", names.arg=c(0:(length(ln)-1)), legend.text=rev(c(\"seed\",\"all\")), col=rev(heat.colors(2)), xlab=\"Linker Degree\", ylab=\"Number of nodes (Log scale)\")\n")
	f.write("axis(2, 0:5, labels=c(0,%s), las=2)\n" % ",".join(map(str, [ 10**i for i in xrange(5) ])))
    else:
	f.write("lA<-matrix(c(ln, lN), nrow=2, byrow=TRUE)\n")
	f.write("barplot(lA, ylim=c(0,round(max(s)/4)), names.arg=c(0:(length(ln)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Linker Degree\", ylab=\"Number of nodes\")\n")

    # Linker degree ratio
    non_seed_ld_ratio_counts = [ 0 ] * 11
    seed_ld_ratio_counts = [ 0 ] * 11
    #for v,ld in node_to_ld.iteritems():
    for v, values in node_to_values.iteritems():
	d, ld, d2, ld2 = values
	#ratio = float(ld) / degrees[v]
	if ld == 0:
	    ratio = 0
	else:
	    ratio = float(ld) / d
	ratio = int(ratio * 10)
	if v in seeds:
	    seed_ld_ratio_counts[ratio] += 1
	else:
	    non_seed_ld_ratio_counts[ratio] += 1
    f.write("rn<-c(%s)\n" % ",".join(map(str, seed_ld_ratio_counts)))
    f.write("rN<-c(%s)\n" % ",".join(map(str, non_seed_ld_ratio_counts)))
    if scale_by_log:
	f.write("rA<-matrix(c(f(rn), f(rN+rn)-f(rn)), nrow=2, byrow=TRUE)\n")
	f.write("barplot(rA, yaxt=\"n\", names.arg=c(0:(length(rn)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Linker degree / Degree\", ylab=\"Number of nodes (Log scale)\")\n")
	f.write("axis(2, 0:5, labels=c(0,%s), las=2)\n" % ",".join(map(str, [ 10**i for i in xrange(5) ])))
    else:
	f.write("rA<-matrix(c(rn, rN), nrow=2, byrow=TRUE)\n")
	f.write("barplot(rA, ylim=c(0,round(max(s)/4)), names.arg=c(0:(length(rn)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Linker degree / Degree\", ylab=\"Number of nodes\")\n")

    # Linker degree ratio at level 2
    non_seed_ld2_ratio_counts = [ 0 ] * 11
    seed_ld2_ratio_counts = [ 0 ] * 11
    for v, values in node_to_values.iteritems():
	d, ld, d2, ld2 = values
	if ld2 == 0:
	    ratio = 0
	else:
	    ratio = float(ld2) / d2
	ratio = int(ratio * 10)
	if v in seeds:
	    seed_ld2_ratio_counts[ratio] += 1
	else:
	    non_seed_ld2_ratio_counts[ratio] += 1
    f.write("rn2<-c(%s)\n" % ",".join(map(str, seed_ld2_ratio_counts)))
    f.write("rN2<-c(%s)\n" % ",".join(map(str, non_seed_ld2_ratio_counts)))
    if scale_by_log:
	f.write("rA2<-matrix(c(f(rn2), f(rN2+rn2)-f(rn2)), nrow=2, byrow=TRUE)\n")
	f.write("barplot(rA2, yaxt=\"n\", names.arg=c(0:(length(rn2)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Linker degree level2 / Degree level2\", ylab=\"Number of nodes (Log scale)\")\n")
	f.write("axis(2, 0:5, labels=c(0,%s), las=2)\n" % ",".join(map(str, [ 10**i for i in xrange(5) ])))
    else:
	f.write("rA2<-matrix(c(rn2, rN2), nrow=2, byrow=TRUE)\n")
	f.write("barplot(rA2, ylim=c(0,round(max(s)/4)), names.arg=c(0:(length(rn2)-1)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Linker degree level2 / Degree level2\", ylab=\"Number of nodes\")\n")

    connected_components = networkx.connected_components(g)
    connected_component_sizes = map(len, connected_components)
    max_size = max(connected_component_sizes)
    connected_component_seed_counts = [ 0 ] * (max_size + 1)
    connected_component_non_seed_counts = [ 0 ] * (max_size + 1)
    for connected_component in connected_components:
	for i in connected_component:
	    if i in seeds:
		connected_component_seed_counts[len(connected_component)] += 1
	    else:
		connected_component_non_seed_counts[len(connected_component)] += 1
    f.write("c<-c(%s)\n" % ",".join(map(str, connected_component_seed_counts)))
    f.write("C<-c(%s)\n" % ",".join(map(str,  connected_component_non_seed_counts)))
    if scale_by_log:
	f.write("cA<-matrix(c(f(c), f(C+c)-f(c)), nrow=2, byrow=TRUE)\n")
	#f.write("barplot(cA, yaxt=\"n\", log=\"x\", xaxp=c(0,length(c),1), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Connected component size (Log scale)\", ylab=\"Number of nodes (Log scale)\")\n")
	f.write("barplot(cA, yaxt=\"n\", xaxt=\"n\", legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Connected component size (Log scale)\", ylab=\"Number of nodes (Log scale)\")\n")
	f.write("axis(2, 0:5, labels=c(0,%s), las=2)\n" % ",".join(map(str, [ 10**i for i in xrange(5) ])))
	f.write("axis(1, seq(1,length(c),round(length(c)/5)), labels=seq(1,length(c),round(length(c)/5)))\n")
    else:
	f.write("cA<-matrix(c(c, C), nrow=2, byrow=TRUE)\n")
	f.write("barplot(cA, ylim=c(0,round(max(s)/4)), names.arg=c(1:length(c)), legend.text=c(\"seed\",\"all\"), col=heat.colors(2), xlab=\"Connected component size\", ylab=\"Number of nodes\")\n")

    f.write("mtext(\'%s\', outer=TRUE, line=-1)\n" % title) 
    f.write("dev.off()\n")
    f.close()
    return


def create_ARFF_network_metrics_file(g, node_to_score, seeds, arff_file_name, calculate_topological_values = False):
    delim = ","
    header = "@RELATION aneurysm\n@ATTRIBUTE id STRING\n@ATTRIBUTE score NUMERIC\n" + \
	        "@ATTRIBUTE degree INTEGER\n@ATTRIBUTE linker_degree INTEGER\n" + \
		"@ATTRIBUTE ld_ratio NUMERIC\n@ATTRIBUTE clustering_coefficient NUMERIC\n" + \
		"@ATTRIBUTE betweenness_centrality NUMERIC\n" + \
	        "@ATTRIBUTE degree2 INTEGER\n@ATTRIBUTE linker_degree2 INTEGER\n" + \
		"@ATTRIBUTE ld_ratio2 NUMERIC\n" + \
		"@ATTRIBUTE class {involved,not-involved}\n@DATA\n" 

    seeds = set(seeds)
    
    if calculate_topological_values:
	print "Calculating betweenness centrality.."
	mapB = networkx.betweenness_centrality(g) 
	##mapB = dict(zip(g.nodes(), range(len(g.nodes())))) 

    if calculate_topological_values:
	print "Calculating clustering coefficients.."
	mapC = networkx.clustering(g, with_labels=True) 

    #print "connected component sizes: ", map(len, networkx.connected_components(g))
    #cliques = networkx.find_cliques(g) # high computational cost
    
    node_to_values = get_node_degree_related_values(g, seeds)

    f = open(arff_file_name, 'w')
    f.write(header)
    for v in g.nodes_iter():
	d, ld, d2, ld2 = node_to_values[v]
	if calculate_topological_values:
	    cc = mapC[v]
	    bc = mapB[v]
	else:
	    cc = 0.0
	    bc = 0.0
	if d == 0: r1 = 0.0
	else: r1 = float(ld)/d
	if d2 == 0: r2 = 0.0
	else: r2 = float(ld2)/d2
	if v in seeds:
	    if node_to_score is not None:
		s=node_to_score[v]
	    else:
		s = 1
	    c="involved"
	else:
	    s="?"
	    c="not-involved"
	# id score degree linker_degree ld_ratio clustering_coeff betweenness_cent d2 ld2 ld_ratio2 class
	# v s d ld n1 cc bc d2 ld2 n2 c
	f.write( ("%s" % delim).join( map(str, [v, s, d, ld, r1, cc, bc, d2, ld2, r2, c]) ) + "\n" )
    f.close()
    return

def create_dot_network_file(g, output_file, seeds=set(), node_to_desc = dict(), ups = set(), downs = set(), draw_type = "all"):
    """
    roots seeds
    remove seed-seed interactions
    colors ups and downs
    """
    #g = create_network_from_sif_file(network_file)
    f = open(output_file, 'w')
    f.write("graph %s {\nforcelabels=false;\noutputorder=edgesfirst\n" % "converted")
    #f.write("NA [style=invis];\n")
    #seeds = set(seed_to_desc.keys()) #set([ line.strip() for line in open(seed_file) ])
    
    linkers = set()
    for node in g.nodes():
	common = set(g.neighbors(node))&seeds
	if len(common) > 1: 
	    linkers.add(node)
	    for i, seed1 in enumerate(common):
		for j, seed2 in enumerate(common):
		    if i<j:
			print node_to_desc[seed1], node_to_desc[node], node_to_desc[seed2]
    #print linkers

    if draw_type == "linker annotated":
	for node in g.nodes():
	    if len(set(g.neighbors(node))&seeds) > 1:
		if node in seeds:
		    f.write("%s [label=\"\" fillcolor=yellow root=true color=yellow height=0.05 width=0.05 shape=rect];\n" % (node)) 
		else:
		    f.write("%s [label=\"\" color=green height=0.05 width=0.05 shape=rect];\n" % (node)) 
	    else:
		if node in seeds:
		    f.write("%s [label=\"\" fillcolor=yellow root=true color=yellow height=0.05 width=0.05 shape=rect];\n" % (node)) 
		else:
		    f.write("%s [label=\"\" fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	ignored = set()
    elif draw_type == "linker_only":
	for node in g.nodes():
	    if node in seeds:
		f.write("%s [label=\"%s\" style=filled fillcolor=yellow root=true color=yellow height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in linkers: 
		if node in ups: 
		    f.write("%s [label=\"\" style=filled fillcolor=green color=green fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
		elif node in downs:
		    f.write("%s [label=\"\" style=filled fillcolor=red color=red fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
		else:
		    f.write("%s [label=\"\" fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	ignored = ((set(g.nodes()) - seeds) - linkers) 
    elif draw_type == "regulated_only":
	for node in g.nodes():
	    if node in seeds:
		f.write("%s [label=\"%s\" style=filled fillcolor=yellow root=true color=yellow height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in ups: 
		f.write("%s [label=\"%s\" style=filled fillcolor=green color=green fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in downs:
		f.write("%s [label=\"%s\" style=filled fillcolor=red color=red fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	ignored = ((set(g.nodes()) - seeds) - ups) - downs 
    elif draw_type == "all":
	for node in g.nodes():
	    if node in seeds:
		f.write("%s [label=\"%s\" style=filled fillcolor=yellow root=true color=yellow height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in ups: 
		f.write("%s [label=\"\" style=filled fillcolor=green color=green fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	    elif node in downs:
		f.write("%s [label=\"\" style=filled fillcolor=red color=red fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	    else:
		f.write("%s [label=\"\" fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	ignored = set()

    for u,v in g.edges_iter():
	if u in seeds and v in seeds:
	    continue
	if u in ignored or v in ignored:
	    continue
	f.write("%s -- %s [color=blue];\n" % (u, v))
    f.write("}\n")
    f.close()
    #os.system("dot %s -Tgif > %s.gif" % (dot_fname, dot_fname))
    return


if __name__ == "__main__":
    main()

    """
    test_network = networkx.Graph()
    test_network.add_nodes_from([1,2,3,4,5,6,7,8,9])
    test_network.add_edges_from([(1,2),(1,3),(1,5),(2,6),(3,4),(5,6),(4,7),(7,8)]) # (1,4)
    test_network2.add_edges_from([(1,2),(1,3),(2,3),(2,4),(4,5)])

    print "original network:"
    print test_network.edges()

    print "preserve topology:"
    random_network = randomize_graph(graph=test_network, randomization_type="preserve_topology")
    print random_network.edges()

    print "preserve topology and node degrees:"
    random_network = randomize_graph(graph=test_network, randomization_type="preserve_topology_and_node_degree")
    print random_network.edges()

    print "preserve degree distribution:"
    randomn = randomize_graph(graph=test_network, randomization_type="preserve_degree_distribution")
    print randomn.edges()

    print "preserve degree distribution and node degree:"
    randomn = randomize_graph(graph=test_network, randomization_type="preserve_degree_distribution_and_node_degree")
    print randomn.edges()

    print "creating big graph..."
    test_power_graph = networkx.Graph()
    test_power_graph.add_nodes_from(range(100000))
    test_power_graph.add_edges_from([ (random.randint(0,99999),random.randint(0,99999)) for x in xrange(1000000) ])
    print "randomizing big network by preserve_topology..."
    randomn = randomize_graph(graph=test_power_graph, randomization_type="preserve_topology")
    print "randomizing by preserve_topology_and_node_degree..."
    randomn = randomize_graph(graph=test_power_graph, randomization_type="preserve_topology_and_node_degree")   
    print "randomizing by preserve_degree_distribution..."
    randomn = randomize_graph(graph=test_power_graph, randomization_type="preserve_degree_distribution")
    """


