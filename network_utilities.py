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

import networkx, random, copy
import os, cPickle, numpy
import types

try:
    from scipy.stats import rankdata
except:
    print "scipy is not installed, rank-based distance methods wont work"

try:
    from external.genrev import NWSteiner, Heuristic, kWalk
except:
    print "GenRev not found, steiner wont work"


#MIN_NUMBER_OF_PERTURBATION = 25
MAX_NUMBER_OF_TRIAL = 10
FINITE_INFINITY = 999999

def main():
    g = create_network_from_sif_file("interactions.sif")
    degrees = dict(g.degree())
    node_to_values = get_node_degree_related_values(g, ["v2","v3"])
    for v in g.nodes():
	print v, degrees[v], node_to_values[v]
    create_R_analyze_network_script(g, seeds = ["v2","v3"], out_path="./")
    g2=prune_graph_at_given_percentage(g,10)
    return

def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
	g = networkx.DiGraph()
    else:
	g = networkx.Graph()
    return g

def create_graph_with_same_type(G):
    return networkx.create_empty_copy(G)

def create_random_graphs(G, n_random, randomization_type, allow_self_edges, out_prefix):
    """
    randomization_type either preserve degree / topology / random
    or prune_k% / permute_k%
    """
    for i in xrange(n_random):
	if randomization_type.startswith("prune_"):
	    percentage = int(randomization_type.split("_")[1])
	    g = prune_graph_at_given_percentage(G, percentage)
	elif randomization_type.startswith("permute_"):
	    percentage = int(randomization_type.split("_")[1])
	    g = permute_graph_at_given_percentage(G, percentage, allow_self_edges = allow_self_edges)
	else:
	    g = randomize_graph(G, randomization_type, allow_self_edges)
	dump_file = "%s_%d.pcl" % (out_prefix, i)
	cPickle.dump(g, open(dump_file, 'w'))
	dump_file = "%s_sp_%d.pcl" % (out_prefix, i)
	sp = get_shortest_path_lengths(g, dump_file)
    return

def get_random_graphs(n_random, out_prefix):
    for i in xrange(n_random):
	dump_file = "%s_%d.pcl" % (out_prefix, i)
	g = cPickle.load(open(dump_file))
	dump_file = "%s_sp_%d.pcl" % (out_prefix, i)
	sp = cPickle.load(open(dump_file))
	yield i, g, sp

def get_number_of_distinct_edges(G):
    edge_list = G.edges()
    edge_set = set()
    for id1, id2, data in edge_list:
        edge_set.add((id1, id2))
    return len(edge_set)

# Imposes restrictions on func, func(a,b) should be always called as func(a, b, None) or func(a, b, dump_file=None)
def dumper(func):
    def wrapper(*args, **kwargs):
	import os
	def modify_kwargs():
	    # remove dump_file argument
	    from copy import deepcopy
	    kwargs_mod = deepcopy(kwargs) 
	    del kwargs_mod["dump_file"]
	    return kwargs_mod
	#print args, kwargs
	if "dump_file" in kwargs:
	    #raise ValueError("dump_file argument not provided")
	    dump_file = kwargs["dump_file"]
	    #kwargs_mod = modify_kwargs()
	    kwargs_mod = kwargs
	    args_mod = args
	else: # assuming the last argument as dump_file
	    dump_file = args[-1]
	    #args_mod = args[:-1]
	    kwargs_mod = kwargs
	    args_mod = args
	#print args_mod, kwargs_mod
	if dump_file is not None:
	    if os.path.exists(dump_file):
	       val = cPickle.load(open(dump_file)) 
	       #print "loading", dump_file
	    else:
		# remove dump_file argument
		#print "dumping", dump_file, kwargs_mod
		val = func(*args_mod, **kwargs_mod)
                if isinstance(val, types.GeneratorType): # newer versions of networkx returns iterators
                    val_new = {source: target_dict for source, target_dict in val}
                    val = val_new
		cPickle.dump(val, open(dump_file, 'w'))
	else:
	    #print "running", kwargs_mod
	    val = func(*args_mod, **kwargs_mod)
	return val
    return wrapper

@dumper
def test(a, b, dump_file=None):
    print a + b
    return

@dumper
def get_shortest_paths(G, dump_file):
    return networkx.shortest_path(G)
    #if dump_file is not None:
    #	if os.path.exists(dump_file):
    #	   sp = cPickle.load(open(dump_file)) 
    #	else:
    #	    sp = networkx.shortest_path(G)
    #	    cPickle.dump(sp, open(dump_file, 'w'))
    #else:
    #	sp = networkx.shortest_path(G)
    #return sp

@dumper
def get_shortest_path_lengths(G, dump_file):
    return networkx.shortest_path_length(G)

@dumper
def get_shortest_path_lengths_subset(G, nodes_subset, dump_file):
    d = networkx.shortest_path_length(G)
    if nodes_subset is None:
	return d
    d_new = dict((node, {}) for node in nodes_subset)
    for node in nodes_subset:
	if node not in d:
	    continue
	d_new[node] = dict((node2, d[node][node2]) for node2 in nodes_subset)
    return d_new

def get_shortest_path_between(G, source_id, target_id):
    return networkx.shortest_path(G, source_id, target_id)

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)

def get_all_shortest_paths_between(G, source_id, target_id):
    return networkx.all_shortest_paths(G, source_id, target_id)

def get_all_paths_from(G, source_id): 
    """
        get all paths from source node to all possible nodes in a dictionary
    """
    return networkx.single_source_dijkstra_path(G, source_id)

def get_all_paths_between(G, source_id, target_id, cutoff=None): 
    return networkx.all_simple_paths(G, source_id, target_id, cutoff)

def get_steiner_tree(G, terminals, sp=None):
    subgraph = NWSteiner.NWConnSteiner(G, terminals, shortestPath=sp)
    return subgraph

def get_heuristic_tree(G, terminals, distance=2, score_increment=0.2):
    # d: search radius, r: the expansion factor
    subgraph = Heuristic.listQuery(G, terminals, d=distance, r=score_increment) # d=1, r=0.2, scoreFun=Heuristic.sum_score
    return subgraph

def get_kwalk_tree(G, terminals):
    subgraph = kWalk.limkWalks(terminals, G) #, L=50, iteration=1)
    return subgraph

@dumper
def get_clustering_coefficient(g, dump_file):
    return networkx.clustering(g)


@dumper
def get_closeness_vitality(g, dump_file):
    return networkx.closeness_vitality(g)


def get_network_radius(g):
    return networkx.radius(g)


def get_network_degree_histogram(g):
    return networkx.degree_histogram(g)


@dumper
def get_node_betweenness(G, dump_file):
    return networkx.betweenness_centrality(G)

@dumper
def get_edge_betweenness(G, dump_file):
    return networkx.edge_betweenness_centrality(G)

def get_edge_betweenness_within_subset(G, subset, edges, consider_alternative_paths=False):
    """
	Calculate edge (relative) betweenness within the paths of a given subset of nodes
	edges is a list of (source, target) pairs for which betweenness is calculated
	Assumes connected component (there is a bath between any two nodes)
	consider_alternative_paths: In case there are multiple shortest paths, the edge is counted if it is in any of them
	Returns a dictionary where values are lists of existance of edge in the shortest path(s) of each edge pair
    """
    edge_to_values = {}
    if consider_alternative_paths:
	for i, u in enumerate(subset):
	    #print i, u, len(subset)
	    for j, v in enumerate(subset):
		if j>=i:
		    continue
		paths_gen = networkx.all_shortest_paths(G, u, v)
		paths = [ [ p for p in path_gen ] for path_gen in paths_gen ]
		for s, t in edges:
		    count = 0
		    total = 0
		    for path in paths:
			#path = [ p for p in path_gen ]
			for k in range(1, len(path)):
			    if (path[k-1], path[k]) == (s, t) or (path[k-1], path[k]) == (t, s):
				count += 1
				break
			total += 1
		    edge_to_values.setdefault((s,t), []).append(float(count)/total)
    else:
	for i, u in enumerate(subset):
	    paths = networkx.shortest_path(G, u)
	    for j, v in enumerate(subset):
		if j>=i:
		    continue
		for s, t in edges:
		    count = 0
		    path = paths[v]
		    for k in range(1, len(path)):
			if (path[k-1], path[k]) == (s, t) or (path[k-1], path[k]) == (t, s):
			    count += 1
			    break
		    edge_to_values.setdefault((s,t), []).append(count)
    edge_to_value = {}
    for edge, values in edge_to_values.iteritems():
	edge_to_value[edge] = float(sum(values)) / len(values)
    return edge_to_value


def get_background_normalized_separation_distance(network, sp, targets, seeds, distance="shortest", n_random=100, n_node_in_bins=100, seed_to_equivalent_nodes = None, seed=None):
    if seed is not None:
	random.seed(seed)
    # Get equivalent nodes from bins
    if seed_to_equivalent_nodes is None:
	bins = get_degree_binning(network, n_node_in_bins)
	seed_to_equivalent_nodes = get_degree_equivalents(seeds, bins, network)
    # Get original score
    d = get_separation(network, sp, targets, seeds, distance)
    values = numpy.empty(n_random)
    for i in xrange(n_random):
	nodes = []
	for seed_node, equivalents in seed_to_equivalent_nodes.iteritems():
	    nodes.append(random.choice(equivalents))
	d_random = get_separation(network, sp, targets, nodes, distance = distance)
	values[i] = d_random
    m, s = numpy.mean(values), numpy.std(values)
    z = d - m
    if s == 0:
	z = 0.0
    else:
	z /= s
    return d, z, (m, s)


def get_source_to_average_target_distance_random(sp, targets, nodes, distance, target_mean_and_std, seed_to_equivalent_nodes, n_random = 100):
    target_to_values = {}
    # Calculate distances of N "equivalent" seeds repeating some hundred times
    for i in xrange(n_random):
	nodes = []
	for seed, equivalents in seed_to_equivalent_nodes.iteritems():
	    nodes.append(random.choice(equivalents))
	target_to_distance_random = get_source_to_average_target_distance(sp, targets, nodes, distance = distance, target_mean_and_std = target_mean_and_std) 
	# Assign the average of all seeds to each target
	for target in targets:
	    target_to_values.setdefault(target, []).append(target_to_distance_random[target])
    return target_to_values

def get_normalized_source_to_target_distance(target_to_distance, target_to_values = None):
    target_to_normalized_distance = {}
    targets = target_to_distance.keys()
    #values = numpy.empty(len(targets))
    for i, target in enumerate(targets):
	if target_to_values is not None: # Normalize using taget_to_values
	    arr = numpy.array(target_to_values[target])
	    m, s = numpy.mean(arr), numpy.std(arr)
	    val = target_to_distance[target] - m
	    if val == 0:
		z = 0.0
	    else:
		z = val / s
	else:
	    z = target_to_distance[target] 
	#values[i] = z
	target_to_normalized_distance[target] = z
    #arr = numpy.sort(values)
    return target_to_normalized_distance #numpy.mean(arr), numpy.std(arr) # numpy.min(arr)


def get_distances_to_node(sp, node):
    if isinstance(sp, dict):
	lengths = sp[node]
    else: # sp is networkx graph
	lengths = networkx.single_source_shortest_path_length(sp, node)
    return lengths


def get_closest_nodes_from_other_set(sp, geneids_source, geneids_target, n_closest=1):
    """
    n_closest: None - get all nodes with min distance
    """
    source_to_target_distance = {}
    for geneid in geneids_source:
	lengths = get_distances_to_node(sp, geneid)
        values = []
        for geneid_target in geneids_target:
            if geneid_target not in lengths:
                #print "Warning: node not connected", geneid_target
                val = FINITE_INFINITY
            else:
                val = lengths[geneid_target]
            values.append((val, geneid_target))
        values.sort()
	if n_closest is None:
	    prev_val = values[0][0]
	    for val, geneid_target in values:
		if val == prev_val:
		    source_to_target_distance.setdefault(geneid, []).append((val, geneid_target))
		else:
		    break
		prev_val = val
	else:
	    source_to_target_distance[geneid] = values[:n_closest]
    return source_to_target_distance


def get_source_to_target_average_distance(network, geneids_source, geneids_target, distance="shortest"):
    """
    Returns average/min distance to target nodes for each source node
    """
    if distance == "tsesolc":
	return get_source_to_target_average_distance(network, geneids_target, geneids_source, distance="closest")
    vals = []
    for geneid in geneids_source:
	values = []
	for geneid_target in geneids_target:
	    val = get_shortest_path_length_between(network, geneid, geneid_target)
	    values.append(val)
	if distance == "shortest":
	    val = numpy.mean(values) 
	elif distance == "closest":
	    val = min(values)
	else:
	    raise ValueError("Unknown distance type " + distance)
	vals.append(val)
    return numpy.mean(vals) 

 
def get_source_to_average_target_distance(sp, geneids_source, geneids_target, distance="shortest", parameters={}, target_mean_and_std = None, exclude_self = False):
    """
    Returns average/min distance to target nodes for each source node
    If target_mean_and_std is provided calculatesn z-score using mean & std of the distance to target from every node
    in addition to the raw distance (first z-score, then raw score)
    """
    if exclude_self: assert geneids_source == geneids_target
    source_to_target_distance = {}
    if distance.startswith("mahalanobis") or distance == "center":
	# find center of seed subnetwork (g-spot)
	center, center_values = get_center_of_subnetwork(sp, geneids_target)
	#center_values = numpy.array(center_values)
	if distance == "mahalanobis-closest" or distance == "mahalanobis-shortest":
	    center_d = numpy.mean(center_values)
	elif distance == "mahalanobis-kernel": 
	    center_d = -numpy.log(numpy.mean([numpy.exp(-value) for value in center_values]))
    elif distance == "jorg.individual":
	target_to_distance = get_source_to_average_target_distance(sp, geneids_target, geneids_target, distance = "closest", parameters = parameters, target_mean_and_std = None, exclude_self=True)
	center_d = numpy.mean(target_to_distance.values())
    elif distance.startswith("rank-"):
	source_to_values = {}
	for geneid in geneids_target:
	    lengths = sp[geneid]
	    ranks = rankdata(lengths.values())
	    for i, geneid_source in enumerate(lengths.keys()):
		if geneid_source in geneids_source:
		    source_to_values.setdefault(geneid_source, []).append(ranks[i])
	for geneid_source, values in source_to_values.iteritems():
	    if distance[len("rank-"):] in ("closest", "closest-min"):
		source_to_target_distance[geneid_source] = numpy.min(values)
	    elif distance[len("rank-"):] == ("shortest", "shortest-min"):
		source_to_target_distance[geneid_source] = numpy.mean(values)
	return source_to_target_distance
    # subsetting seed/target set, e.g. consider only closest 3 seeds/targets with d <= 2 and k >= 10
    if "subset" in parameters:
	geneids_source_mod, geneids_target_mod = get_optimal_subsets_using_parameters(parameters, sp, geneids_source, geneids_target)
	if len(geneids_target_mod) == 0 or len(geneids_source_mod) == 0: 
	    for geneid in geneids_source: 
		val = FINITE_INFINITY
		source_to_target_distance[geneid] = val
	    return source_to_target_distance
        geneids_source = geneids_source_mod
        geneids_target = geneids_target_mod
    for geneid in geneids_source:
	if distance.startswith("net"): # GUILD scores
	    if target_mean_and_std is not None:
		raise ValueError("Normalization is not available for this metric!")
	    val = sp[geneid]
	elif distance == "dsd" or distance.startswith("communicability"):
	    if target_mean_and_std is not None:
		raise ValueError("Normalization is not available for this metric!")
	    DSD, name_to_idx = sp
	    i = name_to_idx[geneid]
	    values = []
	    for geneid_target in geneids_target:
		if exclude_self:
		    if geneid == geneid_target:
			continue
		j = name_to_idx[geneid_target]
		if distance.startswith("communicability"):
		    val = 1 - DSD[i, j]
		else:
		    val = DSD[i, j]
		values.append(val)
	    val = numpy.mean(values) 
	elif distance == "center":
	    val = 0.0
	    for c in center:
		val += sp[geneid][c] # center
	    val /= len(center)
	else:
	    lengths = sp[geneid]
	    values = []
	    for geneid_target in geneids_target:
		if exclude_self:
		    if geneid == geneid_target:
			continue
		if target_mean_and_std is not None:
		    m, s = target_mean_and_std[geneid] # source are drug targets - before [geneid_target] 
		    val = lengths[geneid_target] - m
		    if val == 0:
			val = val
		    else:
			val = val / s  
		else:
		    if geneid_target not in lengths:
			#print "Warning: node not connected", geneid_target
			val = FINITE_INFINITY
		    else:
			val = lengths[geneid_target]
		values.append(val)
	    if exclude_self and len(values) == 0: # for jorg separation metric with sets of one gene
	    	values = [0]
	    if distance == "shortest" or distance == "shortest-min":
		val = numpy.mean(values) 
	    elif distance == "shortest2":
		values = numpy.array(values)
		val = numpy.sqrt(numpy.mean(values*values))
	    elif distance == "kernel" or distance == "kernel-min":
		val = -numpy.log(numpy.mean([numpy.exp(-value-1) for value in values])) 
            elif distance == "kernel2": 
		#val = 0.1/numpy.mean([numpy.exp(-value) for value in values])
		val = numpy.exp(numpy.mean(values))
	    elif distance == "closest" or distance == "closest-min":
		val = min(values)
	    elif distance == "binary":
                val = sum(map(lambda x: x>1, values))
	    elif distance == "closest-adjust":
		values.sort()
		val_min = values[0]
		k = -1
		for val in values:
		    if val == val_min:
			k += 1
		    else:
			break
		val = val_min - k / 100.0 # scale by number of target nodes at that distance
	    elif distance.startswith("knn-"):
		values.sort()
		k = int(distance[-1])
		val = numpy.mean(values[:k])
	    elif distance == "jorg.individual":
		val = min(values) - center_d
	    elif distance == "mahalanobis-closest":
		d = numpy.abs(numpy.array(values) - center_values)
		# val = d # "mahalanobis-x"
		val = numpy.min(d) # "mahalanobis-min"
		# no std scaling "mahalanobis-nostd"
		#val = numpy.mean(d * d) # mahalanobis2
		var = numpy.var(center_values)
		if var > 0:
		#    val /= var
		    val /= numpy.sqrt(var)
		val = numpy.sqrt(val) 
	    elif distance == "mahalanobis-shortest":
		d = numpy.abs(numpy.array(values) - center_values)
		val = numpy.mean(d)
		std = numpy.std(center_values)
		if std == 0:
		    std = 0.01
		val /= std
	    elif distance == "mahalanobis-kernel":
		values = numpy.abs(numpy.array(values) - center_values)
		std = numpy.std(center_values)
		if std == 0:
		    std = 0.01
		values /= std
		val = -numpy.log(numpy.mean([numpy.exp(-value) for value in values]))
	    else:
		raise ValueError("Unknown distance type " + distance)
	source_to_target_distance[geneid] = val
    if distance.startswith("mahalanobis"):
	return source_to_target_distance, center_d
    return source_to_target_distance


def get_optimal_subsets_using_parameters(parameters, sp, geneids_source, geneids_target):
    optimal_subset = parameters["subset"] # "targets" "seeds"
    n_cutoff = FINITE_INFINITY
    m_cutoff = FINITE_INFINITY
    d_cutoff = FINITE_INFINITY
    k_cutoff = 0
    network = None
    geneids_source_sub, geneids_target_sub = geneids_source, geneids_target
    if "n" in parameters:
	n_cutoff = parameters["n"] 
    if "m" in parameters: # used only in both to subset targets
	m_cutoff = parameters["m"] 
    if "d" in parameters:
	d_cutoff = parameters["d"]
    if "k" in parameters:
	k_cutoff = parameters["k"]
	network = parameters["network"]
    if optimal_subset == "seeds": # subseting seeds (geneids_target)
	geneids_target_sub = get_optimal_subset(sp, geneids_source, geneids_target, n_cutoff, d_cutoff, k_cutoff, network) 
    elif optimal_subset == "targets": # subsetting targets (geneids_source)
	geneids_source_sub = get_optimal_subset(sp, geneids_target, geneids_source, n_cutoff, d_cutoff, k_cutoff, network)
    elif optimal_subset == "both": # subsetting both seeds (geneids_target) and targets (geneids_source)
	# subset seeds
	geneids_target_sub = get_optimal_subset(sp, geneids_source, geneids_target, n_cutoff, d_cutoff, k_cutoff, network) 
	# subset targets
	geneids_source_sub = get_optimal_subset(sp, geneids_target, geneids_source, m_cutoff, d_cutoff, k_cutoff, network)
    else:
	raise ValueError("Unknown subsetting type: " + optimal_subset)
    return geneids_source_sub, geneids_target_sub


def get_optimal_subset(sp, geneids_source, geneids_target, n_cutoff=FINITE_INFINITY, d_cutoff=FINITE_INFINITY, k_cutoff=0, network=None):
    """
    Get nodes of geneids_target that are closest (w.r.t. n/d/k restriction) to the geneids_source
    """
    geneids_target_sub = set() 
    values = []
    for geneid_target in geneids_target:
	lengths = get_distances_to_node(sp, geneid_target)
	for geneid in geneids_source:
	    val = lengths[geneid]
	    values.append((val, geneid_target))
    values.sort()
    for val, geneid in values:
	if len(geneids_target_sub) == n_cutoff:
	    break
	if val <= d_cutoff:
	    if network is None:
		geneids_target_sub.add(geneid)
	    else:
		if network.degree(geneid) >= k_cutoff:
		    geneids_target_sub.add(geneid)
    return geneids_target_sub


def get_source_to_average_target_overlap(network, geneids_source, geneids_target, distance):
    """
    Get TOM & MTOM
    """
    source_to_target_distance = {}
    if distance == "tom":
	for geneid in geneids_source:
	    neighbors1 = network.neighbors(geneid)
	    values = []
	    for geneid_target in geneids_target:
		neighbors2 = network.neighbors(geneid_target)
		val = float(len(set(neighbors1)&set(neighbors2)))
		if network.has_edge(geneid, geneid_target):
		    val += 1
		val /= min(len(neighbors1),len(neighbors2)) + 1
		values.append(val)
	    source_to_target_distance[geneid] = 1-numpy.mean(values)
    elif distance == "sttom":
	neighbors1 = set()
	neighbors2 = set()
	for geneid in geneids_source:
	    neighbors1 |= set(network.neighbors(geneid))
	for geneid_target in geneids_target:
	    neighbors2 |= set(network.neighbors(geneid_target))
	val = float(len(set(neighbors1)&set(neighbors2)))
	val /= min(len(neighbors1),len(neighbors2))
	source_to_target_distance["-1"] = 1-val
    elif distance == "mtom":
	count = 0
	for i, geneid_target1 in enumerate(geneids_target):
	    for j, geneid_target2 in enumerate(geneids_target):
		if i < j:
		    if network.has_edge(geneid_target1, geneid_target2):
			count += 1
	for geneid in geneids_source:
	    divisor = None
	    genes = set(geneids_target)
	    genes.add(geneid)
	    for gene1 in genes:
		common_neighbors = set()
		for gene2 in genes:
		    if gene1 == gene2:
			continue
		    neighbors = network.neighbors(gene2)
		    common_neighbors &= set(neighbors)
		if divisor is None or len(common_neighbors) < divisor:
		    divisor = len(common_neighbors)
	    common_neighbors = set(network.neighbors(geneid))
	    for geneid_target in geneids_target:
		neighbors = network.neighbors(geneid_target)
		common_neighbors &= set(neighbors)
	    val = float(len(common_neighbors))
	    for geneid_target in geneids_target:
		if network.has_edge(geneid, geneid_target):
		    val += 1
	    #print val, count, divisor, k
	    val += count
	    k = len(genes)
	    if k == 1:
		val = 1
	    else:
		val /= divisor + (k * (k-1) / 2)
	    source_to_target_distance[geneid] = 1-val 
    elif distance == "mmtom":
	geneids = geneids_source | geneids_target
	count = 0
	for i, geneid1 in enumerate(geneids):
	    for j, geneid2 in enumerate(geneids):
		if i < j:
		    if network.has_edge(geneid1, geneid2):
			count += 1
	divisor = None 
	common_neighbors = set()
	for geneid1 in geneids:
	    common_neighbors_inner = set()
	    for geneid2 in geneids:
		if geneid1 == geneid2:
		    continue
		neighbors = network.neighbors(geneid2)
		if len(common_neighbors_inner) == 0:
		    common_neighbors_inner = set(neighbors)
		else:
		    common_neighbors_inner &= set(neighbors)
	    if divisor is None or len(common_neighbors_inner) < divisor:
		divisor = len(common_neighbors_inner)
	    neighbors = network.neighbors(geneid1)
	    if len(common_neighbors) == 0:
		common_neighbors = set(neighbors)
	    else:
		common_neighbors &= set(neighbors)
	val = float(len(common_neighbors))
	#print val, count, divisor, k
	val += count
	k = len(geneids)
	val /= divisor + (k * (k-1) / 2)
	source_to_target_distance["-1"] = 1-val 
    else:
	raise ValueError("Unknown distance type " + distance)
    return source_to_target_distance


def get_separation(network, sp, targets, seeds, distance, parameters={}, averaging_function=lambda x: numpy.mean(x)): 
    """
    Potential averaging functions
    val = numpy.min(x) 
    val = numpy.mean(x) 
    val = -numpy.log(numpy.mean([numpy.exp(-value-1) for value in x]))
    jorg / mmtom / avg distances: knn-x, shortest, kernel, mahalanobis, tom, mtom, ... 
    """
    if distance.startswith("jorg-"):
	distance = distance[len("jorg-"):] # closest was default before
	target_to_distance = get_source_to_average_target_distance(sp, seeds, seeds, distance = distance, parameters = parameters, target_mean_and_std = None, exclude_self=True)
	values = target_to_distance.values()
	d1 = averaging_function(values)
	target_to_distance = get_source_to_average_target_distance(sp, targets, targets, distance = distance, parameters = parameters, target_mean_and_std = None, exclude_self=True)
	values = target_to_distance.values()
	d2 = averaging_function(values)
	target_to_distance = get_source_to_average_target_distance(sp, targets, seeds, distance = distance, parameters = parameters, target_mean_and_std = None)
	values = target_to_distance.values()
	target_to_distance = get_source_to_average_target_distance(sp, seeds, targets, distance = distance, parameters = parameters, target_mean_and_std = None)
	values.extend(target_to_distance.values())
	d12 = averaging_function(values)
	val = d12 - (d1 + d2) / 2.0
    elif distance.startswith("mahalanobis-jorg-"):
	distance = distance[len("mahalanobis-jorg-"):] # mahalanobis is default before
	target_to_distance, d1 = get_source_to_average_target_distance(sp, targets, seeds, distance, parameters = parameters)
	values = target_to_distance.values()
	target_to_distance, d2 = get_source_to_average_target_distance(sp, seeds, targets, distance, parameters = parameters)
	values.extend(target_to_distance.values())
	d12 = averaging_function(values)
	val = d12 - (d1 + d2) / 2.0
    elif distance == "mahalanobis-pairwise":
	target_to_distance, center_d = get_source_to_average_target_distance(sp, targets, seeds, distance = "mahalanobis-shortest", parameters = parameters)
	values = target_to_distance.values()
	target_to_distance, center_d = get_source_to_average_target_distance(sp, seeds, targets, distance = "mahalanobis-shortest", parameters = parameters)
	values.extend(target_to_distance.values())
	val = averaging_function(values)
    elif distance == "center-pairwise":
	center_targets, center_values = get_center_of_subnetwork(sp, targets)
	center_seeds, center_values = get_center_of_subnetwork(sp, seeds)
	values = []
	for c_t in center_targets:
	    for c_s in center_seeds:
		values.append(sp[c_t][c_s])
	val = averaging_function(values)
    elif distance == "closest-pairwise":
	values = []
	for geneid in targets:
	    lengths = sp[geneid]
	    inner_values = []
	    for geneid_seed in seeds:
		if geneid == geneid_seed:
		    val = 0
		else:
		    val = lengths[geneid_seed]
		inner_values.append(val)
	    values.append(min(inner_values))
	for geneid in seeds:
	    lengths = sp[geneid]
	    inner_values = []
	    for geneid_target in targets:
		if geneid == geneid_target:
		    val = 0
		else:
		    val = lengths[geneid_target]
		inner_values.append(val)
	    values.append(min(inner_values))
	val = averaging_function(values)
    elif distance == "shortest-pairwise":
	values = []
	for geneid in targets:
	    lengths = sp[geneid]
	    for geneid_seed in seeds:
		if geneid == geneid_seed:
		    val = 0
		else:
		    val = lengths[geneid_seed]
		values.append(val)
	val = averaging_function(values)
    elif distance == "dsd-pairwise" or distance == "communicability-pairwise":
	DSD, name_to_idx = sp
	values = []
	for geneid in targets:
	    i = name_to_idx[geneid]
	    for geneid_seed in seeds:
		if geneid == geneid_seed:
		    val = 0
		else:
		    j = name_to_idx[geneid_seed]
		    if distance.startswith("communicability"):
			val = 1 - DSD[i, j]
		    else:
			val = DSD[i, j]
		values.append(val)
	for geneid in seeds:
	    i = name_to_idx[geneid]
	    for geneid_target in targets:
		if geneid == geneid_target:
		    val = 0
		else:
		    j = name_to_idx[geneid_target]
		    val = DSD[i, j]
		values.append(val)
	val = averaging_function(values)
    elif distance == "kernel-pairwise":
	values = []
	for geneid in targets:
	    lengths = sp[geneid]
	    for geneid_seed in seeds:
		if geneid == geneid_seed:
		    val = 0
		else:
		    val = lengths[geneid_seed]
		values.append(val)
	val = -numpy.log(numpy.sum([numpy.exp(-value-1) for value in values])) / len(values)
    else:
	if distance.endswith("tom") or distance.endswith("tom-min"):
	    target_to_distance = get_source_to_average_target_overlap(network, targets, seeds, distance)
	elif distance.startswith("mahalanobis"):
	    target_to_distance, center_d = get_source_to_average_target_distance(sp, targets, seeds, distance, parameters = parameters)
	elif distance == "tsesolc":
	    target_to_distance = get_source_to_average_target_distance(sp, seeds, targets, "closest", parameters = parameters)
	else:
	    target_to_distance = get_source_to_average_target_distance(sp, targets, seeds, distance, parameters = parameters)
	values = target_to_distance.values()
	if distance.endswith("-min"):
	    val = numpy.min(values)
	else:
	    val = averaging_function(values)
    return val

def get_adjacency_matrix(G):
    nodes = G.nodes() 
    n = len(nodes)
    node_to_idx = {}
    idx_to_node = {}
    for i, node in enumerate(nodes):
	node_to_idx[node] = i
	idx_to_node[i] = node
    adjacency = numpy.zeros((n, n))
    for u, v in G.edges():
	i = node_to_idx[u]
	j = node_to_idx[v]
	if i == j:
	    continue
	adjacency[i, j] = 1
	adjacency[j, i] = 1
    return adjacency, idx_to_node

@dumper
def get_communicability_distances(G, degree_scaling_function, dump_file):
    """
    G - networkx graph object #, assumes graph is fully connected.

    returns D distance matrix 
    """
    from scipy import linalg
    #from scipy.sparse import linalg # need to update for expm 
    adjacency, idx_to_node = get_adjacency_matrix(G)
    n = numpy.size(adjacency[0])
    degree = numpy.zeros((n, 1))
    for i in xrange(0, n):
	degree[i] = sum(adjacency[i])
    adjacency /= degree_scaling_function(degree)
    D = linalg.expm(adjacency)
    return D, idx_to_node 


@dumper
def get_difusion_state_distances(G, nRW, dump_file):
    '''
    modified from
    calcDSD.py -- This module parse calculates DSD given the adjacency matrix
	and output accoding to options

    DSD version 0.5, Copyright (C) 2013, Tufts University
    @author -- Mengfei Cao, mcao01@cs.tufts.edu
    161 College Ave., Medford, MA 02155, USA

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation version 2.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
    MA 02110-1301, USA

    '''
    """
    G - networkx graph object, assumes graph is fully connected.

    nRW - the length of random walks used to calculate DSD
          if nRW = -1, then calculate 

    returns DSD score matrix 
    """
    adjacency, idx_to_node = get_adjacency_matrix(G)
    #### p for transition matrix
    n = numpy.size(adjacency[0])
    p = numpy.zeros((n, n))
    degree = numpy.zeros((n, 1))
    for j in xrange(0, n):
        degree[j] = sum(adjacency[j])
        for i in xrange(0, n):
            if degree[j] != 0:
                p[j] = adjacency[j]/degree[j]

    if nRW >= 0:
        #### c for visit count matrix
        #### for example, c(2,3) is the number of times
        ####     that node 3 is visited via random walks
        ####     starting from node 2
        c = numpy.eye(n)
        for rw in xrange(0, nRW):
            #print c
            c = numpy.dot(c, p) + numpy.eye(n)
            #### this is the c matrix for random walks with
            #### length rw, i.e. if length is 0, then c is
            #### simply the identiy matrix, visiting itself
            #### if length is 1, then c is eye() plus the
            #### the one step transition matrix
    else:
        c = numpy.eye(n)
        c = c - p
        pi = (degree.conj().T)/sum(degree)
        c = c + numpy.tile(pi, (n, 1))
        c = numpy.linalg.inv(c)

    DSD = numpy.zeros((n, n))
    for i in xrange(0, n):
        for j in xrange(i+1, n):
            if degree[i] and degree[j]:
                DSD[i, j] = numpy.linalg.norm((c[i, :]-c[j, :]), ord=1)
                DSD[j, i] = DSD[i, j]
            else:
                DSD[i, j] = -1
                DSD[j, i] = -1
    return DSD, idx_to_node 

# Below not used due to memory problems 
def convert_distance_matrix_to_dictionary(D, idx_to_node):
    sp = {}
    for i in xrange(0, n):
	u = idx_to_node[i]
	sp[u] = {}
        for j in xrange(0, n):
	    val = D[i, j]
            if val < 0:
		continue
	    v = idx_to_node[j]
	    sp[u][v] = val
    return sp 


def get_center_of_subnetwork(sp, nodes):
    center = None
    center_d = None
    center_values = None
    for i, node1 in enumerate(nodes):
	lengths = sp[node1]
	values = []
	for j, node2 in enumerate(nodes):
	    if i == j:
		values.append(0)
	    else:
		values.append(lengths[node2])
	d = sum(values)
	if center_d is None or d < center_d:
	    center = [ node1 ]
	    center_d = d
	    #center_values = [ values ]
	    center_values = numpy.array(values, dtype=float)
	elif d == center_d:
	    center.append(node1)
	    center_values += numpy.array(values)
	    #center_values.append(values)
    center_values /= len(center)
    return center, center_values


def get_degree_binning(g, bin_size, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree(): #.iteritems(): # iterator in networkx 2.0
	if lengths is not None and node not in lengths:
	    continue
	degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    values.sort()
    bins = []
    i = 0
    while i < len(values):
	low = values[i]
	val = degree_to_nodes[values[i]]
	while len(val) < bin_size:
	    i += 1
	    if i == len(values):
		break
	    val.extend(degree_to_nodes[values[i]])
	if i == len(values):
	    i -= 1
	high = values[i]
	i += 1 
	#print i, low, high, len(val) 
	if len(val) < bin_size:
	    low_, high_, val_ = bins[-1]
	    bins[-1] = (low_, high, val_ + val)
	else:
	    bins.append((low, high, val))
    return bins


def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
	d = g.degree(seed)
	for l, h, nodes in bins:
	    if l <= d and h >= d:
		mod_nodes = list(nodes)
		mod_nodes.remove(seed)
		seed_to_nodes[seed] = mod_nodes
		break
    return seed_to_nodes


def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False, seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
	random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in xrange(n_random):
	if degree_aware:
	    if connected:
		raise ValueError("Not implemented!")
	    nodes_random = set()
	    node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
	    for node, equivalent_nodes in node_to_equivalent_nodes.iteritems():
	    	#nodes_random.append(random.choice(equivalent_nodes))
		chosen = random.choice(equivalent_nodes)
		for k in xrange(20): # Try to find a distinct node (at most 20 times)
		    if chosen in nodes_random:
			chosen = random.choice(equivalent_nodes)
	    	nodes_random.add(chosen)
	    nodes_random = list(nodes_random)
	else:
	    if connected:
		nodes_random = [ random.choice(nodes) ]
		k = 1
		while True:
		    if k == len(nodes_selected):
			break
		    node_random = random.choice(nodes_random)
		    node_selected = random.choice(network.neighbors(node_random))
		    if node_selected in nodes_random:
			continue
		    nodes_random.append(node_selected)
		    k += 1
	    else:
		nodes_random = random.sample(nodes, len(nodes_selected))
	values.append(nodes_random)
    return values


def get_pairwise_distances_between_nodes(g, sources, targets=None):
    symetric = False
    if targets is None:
	targets = sources
	symetric = True
    values = []
    for i, u in enumerate(sources):
	sp_lengths = networkx.single_source_dijkstra_path_length(g, u)
	for j, v in enumerate(targets):
	    if symetric and j>=i:
		continue
	    values.append(sp_lengths[v])
    return values


def get_interacting_pairs(g, nodes1, nodes2):
    values = []
    for node1 in nodes1:
	for node2 in nodes2:
	    if node1 == node2 or g.has_edge(node1, node2):
		values.append((node1, node2))
    return values


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
                        new_graph.add_edge(listNodes[x],listNodes[y]) #,len(sp)-1)
            
    return new_graph


def get_subgraph(G, nodes):
    """
	NetworkX subgraph method wrapper
    """
    return G.subgraph(nodes)


def is_subgraph(G, H):
    """
	Checks whether G is a subgraph of H, that is whether all the edges of G belong to H
    """
    edges1 = set(G.edges())
    edges2 = set(H.edges())
    return len(edges1) == len(edges1 & edges2)


def merge_graphs(G, H):
    """
	NetworkX union method wrapper
    """
    return networkx.union(G, H)

def get_neighborhood_subgraph(g, nodes):
    """
	Returns the subgraph induced by the given nodes and their neighbors
    """
    neighbors = []
    for node in nodes:
	neighbors.extend(g.neighbors(node))
	neighbors.extend(nodes)
	g_sub = g.subgraph(neighbors)
    return g_sub

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
	result_list = [c for c in sorted(networkx.connected_components(G), key=len, reverse=True)]

    return result_list


def create_network_from_sif_file(network_file_in_sif, use_edge_data = False, delim = None, include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type = use_edge_data, delim = delim)
    g = create_graph()
    if include_unconnected:
	g.add_nodes_from(setNode)
    if use_edge_data:
	for e,w in dictEdge.iteritems():
	    u,v = e
	    g.add_edge(u,v,w=w) #,{'w':w})
    else:
	g.add_edges_from(setEdge)
    return g


def create_network_from_first_two_columns(network_file, delim = None):
    g = create_graph()
    for line in open(network_file):
	id1, id2 = line.strip().split(delim)[:2]
	g.add_edge(id1, id2)
    return g


def output_node_info(g, node_to_label, out_file, default_label="", node_to_name=None):
    f = open(out_file, 'w')
    f.write("Column 1\tColumn 2\n")
    for u in g.nodes():
	label = default_label
	if u in node_to_label:
	    label = node_to_label[u]
	name = u
	if node_to_name is not None:
	    if u in node_to_name:
		name = node_to_name[u]
	f.write("%s\t%s\n" % (name, label))
    f.close()

def output_network_in_sif(g, output_file_name, node_to_desc=None, delim = " ", include_unconnected=True, remove_self=True):
    f = open(output_file_name, 'w')
    included_nodes = set()
    for u,v in g.edges():
	try:
	    weight = str(g.get_edge_data(u,v)['w'])
	except:
	    weight = "1"
	if node_to_desc is not None:
	    desc1, desc2 = node_to_desc[u], node_to_desc[v]
	else:
	    desc1, desc2 = u, v
	if remove_self:
	    if desc1 == desc2:
		continue
	f.write("%s%s%s%s%s\n" % (desc1, delim, weight, delim, desc2) )
	included_nodes.add(desc1)
	included_nodes.add(desc2)
    if include_unconnected:
	for u in g.nodes():
	    if node_to_desc is not None:
		desc = node_to_desc[u]
	    else:
		desc = u
	    if desc not in included_nodes:
		f.write("%s\n" % desc)
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
    degrees = zip(*list(g.degree()))[1]
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
    degrees = g.degree()
    subgraph_nodes = []
    for u, d in degrees: 
	if degree_threshold is None or d <= degree_threshold:
	    subgraph_nodes.append(u)
    g_filtered = g.subgraph(subgraph_nodes)
    if largest_connected_component:
	component_nodes = max(networkx.connected_components(g_filtered), key=len).nodes()
	g_filtered = g_filtered.subgraph(component_nodes)
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
    flag = False
    f=open(file_name)
    for line in f:
	if delim is None:
	    words = line.rstrip("\n").split()
	else:
	    words = line.rstrip("\n").split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
	    if data_to_float:
		score = float(words[1])
	    else:
		score = words[1]
            dictNode[id1] = score
        elif len(words) >= 3: 
	    if len(words) > 3:
		flag = True
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
    if flag:
	print "Warning: Ignored extra columns in the file!"
    return setNode, setEdge, dictNode, dictEdge


def get_jaccard_index_map(g):
    edge_to_jaccard = {}
    for u,v in g.edges():
	u_neighbors = set(g.neighbors(u))
	v_neighbors = set(g.neighbors(v))
	edge_to_jaccard[(u,v)] = float(len(u_neighbors & v_neighbors)) / len(u_neighbors | v_neighbors)
    return edge_to_jaccard 


def randomize_graph(graph, randomization_type, allow_self_edges = False):
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

    if randomization_type == "same_degree_sequence":
	# Takes ages to find a suitable conformation for large graphs
	sequence = graph.degree().values()
	new_graph = None
	while new_graph is None: 
	    new_graph = networkx.random_degree_sequence_graph(sequence)
	return new_graph

    if randomization_type == "graph_tool_correlated":
	try:
	    import graph_tool
	except:
	    raise ValueError("Graph tool package not installed")
	    return
	new_graph = graph.copy()
	graph_tool.generation.random_rewire(new_graph, model='uncorrelated', n_iter=1, edge_sweep=True, parallel_edges=False, self_loops=False, vertex_corr=None, block_membership=None, alias=True, cache_probs=True, persist=False, ret_fail=False, verbose=False)
	return new_graph

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
		new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
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
	    new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
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
	degree_map = dict(networkx.degree(new_graph))
	nodes = new_graph.nodes()

	# Map graph from random model to new graph
        for edge in new_graph.edges():
	    if len(available_edges) > 0:
		edge_org = available_edges.pop()
		if debug:
		    print "From random:", (edge[0], edge[1])
		new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
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
	    new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
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
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge[0],edge[1]))
        
    elif randomization_type=="preserve_topology": # shuffle_nodes
        nodes = graph.nodes()
        random_nodes = graph.nodes()
        random.shuffle(random_nodes)
        equivalences = dict([(nodes[i],random_nodes[i]) for i in xrange(len(nodes))])
        new_graph.add_edges_from([ (equivalences[current_edge[0]], equivalences[current_edge[1]], graph.get_edge_data(current_edge[0],current_edge[1])) for current_edge in graph.edges() ])

    elif randomization_type=="preserve_topology_and_node_degree": # shuffle_nodes_within_same_degree
        nodes_by_degree = dict( (degree,[]) for u, degree in graph.degree() ) #.values()
        graph_degree = dict(graph.degree())
        [ nodes_by_degree[graph_degree[node]].append(node) for node in graph_degree ]
        equivalences = {}
        for current_degree in nodes_by_degree.keys():
            nodes = nodes_by_degree[current_degree]
            random_nodes = list(nodes)
            random.shuffle(random_nodes)
            equivalences.update(dict([(nodes[i],random_nodes[i]) for i in xrange(len(nodes))]))
        new_graph.add_edges_from([ (equivalences[current_edge[0]],equivalences[current_edge[1]], graph.get_edge_data(current_edge[0],current_edge[1])) for current_edge in graph.edges() ])
        
    elif randomization_type=="preserve_degree_distribution":
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge_data(current_node1, current_node2))
        max_degree = sorted(zip(*list(graph.degree()))[1])[-1] #.values()
        nodes_by_degree = dict( (degree, {}) for degree in xrange(max_degree+1) )
        graph_degree = dict(graph.degree())
        [ nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree ]
        n_perturbation = random.randint(2*n_edge/3, n_edge) # Perturb at least 66% of the edges
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
		if target_id == source_id:
		    target_degree -= 1
                del nodes_by_degree[target_degree][target_id] 
                nodes_by_degree[target_degree-1].setdefault(target_id)
                ## not very important to check for cases where new_source = source (v.v. for targets) 
                new_target_id = random.choice(nodes_by_degree[target_degree-1].keys())
		if source_id == target_id:
		    new_source_id = new_target_id
		else:
		    new_source_id = random.choice(nodes_by_degree[source_degree-1].keys())
		if debug:
		    print source_id, target_id, " / ", new_source_id, new_target_id
		    print source_degree, target_degree
                ## check if going to add an existing edge or self edge
                if new_graph.has_edge(new_source_id, new_target_id) or (not allow_self_edges and new_source_id == new_target_id):
                    del nodes_by_degree[target_degree-1][target_id] 
                    nodes_by_degree[target_degree].setdefault(target_id)
                    del nodes_by_degree[source_degree-1][source_id] 
                    nodes_by_degree[source_degree].setdefault(source_id)
                    continue
		if debug:
		    print "rm %s %s" % (source_id, target_id)
                edge_data = new_graph.get_edge_data(source_id, target_id)
                new_graph.remove_edge(source_id, target_id)
		if debug:
		    print "add %s %s" % (new_source_id, new_target_id)
                new_graph.add_edge(new_source_id, new_target_id, edge_data)
                del nodes_by_degree[target_degree-1][new_target_id] 
                nodes_by_degree[target_degree].setdefault(new_target_id)
		if new_source_id == new_target_id and source_id != target_id:
		    source_degree += 1
                del nodes_by_degree[source_degree-1][new_source_id] 
                nodes_by_degree[source_degree].setdefault(new_source_id)
                break
        randomize_graph(new_graph, "preserve_topology")

    elif randomization_type=="preserve_degree_distribution_and_node_degree":
        ## add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge_data(current_node1, current_node2))
        nodes_by_degree = dict( (degree,{}) for u, degree in graph.degree() )
        graph_degree = dict(graph.degree())
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
                edge_data_1 = new_graph.get_edge_data(source_id, target_id)
                edge_data_2 = new_graph.get_edge_data(new_source_id, new_target_id)
                new_graph.remove_edge(source_id, target_id)
                new_graph.remove_edge(new_source_id, new_target_id)
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
	new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge[0],edge[1]))
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
    #degrees = g.degree() 
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
    for v in g.nodes():
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

def create_dot_network_file(g, output_file, seeds=set(), node_to_desc = dict(), ups = set(), downs = set(), weak_edges = set(), draw_type = "all"):
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

    ignored = set()
    if draw_type == "linker annotated":
	for node in g.nodes():
	    if len(set(g.neighbors(node))&seeds) > 1:
		if node in seeds:
		    f.write("%s [label=\"\" fillcolor=red root=true color=red height=0.05 width=0.05 shape=rect];\n" % (node)) 
		else:
		    f.write("%s [label=\"\" color=green height=0.05 width=0.05 shape=rect];\n" % (node)) 
	    else:
		if node in seeds:
		    f.write("%s [label=\"\" fillcolor=red root=true color=red height=0.05 width=0.05 shape=rect];\n" % (node)) 
		else:
		    f.write("%s [label=\"\" fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	ignored = set()
    elif draw_type == "linker_only":
	for node in g.nodes():
	    if node in seeds:
		f.write("%s [label=\"%s\" style=filled fillcolor=red root=true color=red height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in linkers: 
		if node in ups: 
		    f.write("%s [label=\"\" style=filled fillcolor=green color=green fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
		elif node in downs:
		    f.write("%s [label=\"\" style=filled fillcolor=blue color=blue fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
		else:
		    f.write("%s [label=\"\" fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	ignored = ((set(g.nodes()) - seeds) - linkers) 
    elif draw_type == "regulated_only":
	for node in g.nodes():
	    if node in seeds:
		f.write("%s [label=\"%s\" style=filled fillcolor=red root=true color=red height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in ups: 
		f.write("%s [label=\"%s\" style=filled fillcolor=green color=green fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in downs:
		f.write("%s [label=\"%s\" style=filled fillcolor=blue color=blue fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	ignored = ((set(g.nodes()) - seeds) - ups) - downs  
    elif draw_type == "seeds_only":
	for node in g.nodes():
	    if node in seeds:
		f.write("%s [label=\"%s\" style=filled fillcolor=red root=true color=red height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	ignored = set(g.nodes()) - seeds  
    elif draw_type == "all":
	for node in g.nodes():
	    if node in seeds:
		f.write("%s [label=\"%s\" style=filled fillcolor=red root=true color=red height=0.05 width=0.05 shape=rect];\n" % (node, node_to_desc[node])) 
	    elif node in ups: 
		f.write("%s [label=\"\" style=filled fillcolor=green color=green fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	    elif node in downs:
		f.write("%s [label=\"\" style=filled fillcolor=blue color=blue fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node)) 
	    elif not any(map(lambda x: x!=node and x in node_to_desc, g.neighbors(node))):
		continue # skip nodes with only self edge 
	    else:
		f.write("%s [label=\"%s\" fixedsize=true height=0.05 width=0.05 shape=rect];\n" % (node, "")) # node_to_desc[node]
	ignored = set()

    for u,v in g.edges():
	if u in ignored or v in ignored:
	    continue
	if u == v: # skip self edges
	    continue
	# draw dashed lines for weak edges
	if len(weak_edges) > 0:
	    if (u, v) in weak_edges:
		f.write("%s -- %s [color=darkgrey, style=dashed];\n" % (u, v))
	    else:
		f.write("%s -- %s [color=blue];\n" % (u, v))
	else:
	    # skip edges between seeds
	    #if u in seeds and v in seeds:
	    #	continue
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


