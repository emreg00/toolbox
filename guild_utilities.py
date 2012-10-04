
from toolbox import network_utilities, file_converter, stat_utilities
import os

def main():
    network_file = "/home/emre/arastirma/data/collaboration/billur/9606/network_no_tap_geneid.sif"
    #seed_file = "/home/emre/arastirma/data/collaboration/billur/brain_seeds_geneid.txt"
    seed_file = "/home/emre/arastirma/data/collaboration/billur/lung_seeds_geneid.txt"
    scoring_folder = "./test/"
    executable_dir = "/home/emre/arastirma/netzcore/src/"
    prepare_scoring(network_file, seed_file, scoring_folder, non_seed_score=0.01, seed_score=1.0, edge_score=1.0, n_sample=100, delim=" ")
    #run_scoring(scoring_folder, executable_dir, scoring_type="netzcore", parameters={"n_iteration":5, "n_sample":100, "sampling_prefix":scoring_folder+"sampled_graph."}, qname=None)
    run_scoring(scoring_folder, executable_dir, scoring_type="netcombo")
    return

def prepare_scoring(network_file, seed_file, scoring_folder="./", non_seed_score=0.01, seed_score=1.0, edge_score=1.0, n_sample=100, delim=" "):
    """
	network file: sif-like format where edge type is edge score: A 0.5 B
	seed file: sif-like format where nodes and their score are given: A 0.1
    """
    if not os.path.exists(scoring_folder):
	os.mkdir(scoring_folder)
    # Read node info from network file (use network file as edge file)
    print "Creating edge score file"
    edge_file = scoring_folder + "edge_scores.sif" #network_file.split("/")[-1] + ".converted"
    if os.path.exists(edge_file):
	print "\tEdge file exists, overwriting!"
    nodes, edges, dummy, edge_to_data = network_utilities.get_nodes_and_edges_from_sif_file(network_file, store_edge_type = True, delim = delim, data_to_float=False)
    edge_to_weight = {}
    f = open(edge_file, 'w')
    for edge in edges:
	data = edge_to_data[edge]
	try:
	    score = float(data)
	except:
	    score = edge_score
	edge_to_weight[edge] = score
	f.write("%s%s%f%s%s\n" % (edge[0], delim, score, delim, edge[1]))
    f.close()
    # Create node file (ignore seeds that are not in the network and assign non-seed scores)
    print "Creating node score file"
    from random import shuffle
    node_file = scoring_folder +  "node_scores.sif" #seed_file.split("/")[-1] + ".converted"
    seeds, dummy, seed_to_data, dummy = network_utilities.get_nodes_and_edges_from_sif_file(seed_file, store_edge_type = False, delim = delim, data_to_float=False)
    f = open(node_file, 'w')
    node_to_data = {}
    for node in nodes:
	if node in seeds:
	    if seed_to_data is not None:
		score = seed_to_data[node]
	    else:
		score = seed_score
	else:
	    score = non_seed_score
	node_to_data[node] = score
	f.write("%s%s%f\n" % (node, delim, score))
    f.close()
    # Create background node file (selects k non-seeds randomly where k is the number of seeds)
    print "Creating background node score file"
    non_seeds = list(nodes - seeds)
    shuffle(non_seeds)
    random_seeds = set(non_seeds[:len(seeds)])
    bg_node_file = scoring_folder +  "node_scores_background.sif" #seed_file.split("/")[-1] + ".converted"
    f = open(bg_node_file, 'w')
    if seed_to_data is not None: seed_scores = seed_to_data.values()
    for node in nodes:
	if node in random_seeds:
	    if seed_to_data is not None:
		score = seed_scores.pop()
	    else:
		score = seed_score
	else:
	    score = non_seed_score
	f.write("%s%s%f\n" % (node, delim, score))
    f.close()
    # Create modified edge file using node scores for netshort
    print "Creating node score converted edge file (for netshort)"
    nd_edge_file = scoring_folder + "edge_scores_netshort.sif" #network_file.split("/")[-1] + ".converted_for_netshort"
    f = open(nd_edge_file, 'w')
    for u,v in edges:
	score_u = node_to_data[u]
	score_v = node_to_data[v]
	weight = edge_to_weight[(u, v)]
	f.write("%s%s%f%s%s\n" % (u, delim, weight*(score_u + score_v) / 2, delim, v))
    f.close()
    # Create random network files for netzcore
    print "Creating random networks (for netzcore)"
    sampling_prefix = scoring_folder + "sampled_graph."
    if os.path.exists(sampling_prefix+"%s"%n_sample):
	print "\tSampled networks exists, skipping this step!"
    else:
	g = network_utilities.create_network_from_sif_file(network_file_in_sif = edge_file, use_edge_data = True, delim = delim)
	for i in xrange(1,n_sample+1):
	    g_sampled = network_utilities.randomize_graph(graph=g, randomization_type="preserve_topology_and_node_degree")
	    network_utilities.output_network_in_sif(g_sampled, sampling_prefix+"%s"%i)
    return

def run_scoring(scoring_folder, executable_dir, scoring_type="netscore", parameters={"n_iteration":2, "n_repetition":3, "n_sample":100, "sampling_prefix":"./sampled_graph."}, qname=None):
    """
    scoring_type: netscore | netzcore | netshort | netcombo
    qname: sbi | sbi-short | bigmem
    """

    def score(scoring_type, qname, node_file, edge_file, output_file, parameters):
	output_file += ".%s" % scoring_type
	if scoring_type == "netscore": 
	    score_command = executable_dir + "scoreNetwork/scoreN -s s -n %s -e %s -o %s -r %d -i %d" % (node_file, edge_file, output_file, parameters["n_repetition"], parameters["n_iteration"])
	elif scoring_type == "netzcore": 
	    score_command = executable_dir + "scoreNetwork/scoreN -s z -n %s -e %s -o %s -i %d -x %d -d %s" % (node_file, edge_file, output_file, parameters["n_iteration"], parameters["n_sample"], parameters["sampling_prefix"])
	elif scoring_type == "netshort": 
	    score_command = executable_dir + "scoreNetwork/scoreN -s d -n %s -e %s -o %s" % (node_file, parameters["nd_edge_file"], output_file)
	else:
	    raise ValueError("Invalid scoring type!")
	if qname is None:
	    os.system(score_command)
	else:
	    os.system("qsub -cwd -o out -e err -q %s -N %s -b y %s" % (qname, scoring_type, score_command))
	return

    edge_file = scoring_folder + "edge_scores.sif" 
    node_file = scoring_folder +  "node_scores.sif" 
    bg_node_file = scoring_folder +  "node_scores_background.sif" 
    nd_edge_file = scoring_folder + "edge_scores_netshort.sif"
    sampling_prefix = scoring_folder + "sampled_graph."
    output_file = scoring_folder + "output_scores.sif"
    bg_output_file = scoring_folder + "output_scores_background.sif"
    if not os.path.exists(node_file) or not os.path.exists(edge_file):
	print "Input files not found!\nMake sure that you have run prepare_scoring first and that you provide the correct path."
	return
    # Run scoring algorithm
    parameters["sampling_prefix"] = sampling_prefix
    if scoring_type == "netcombo":
	scoring = "netscore"
	parameters={"n_repetition":3, "n_iteration":2}
	score(scoring, qname, node_file, edge_file, output_file, parameters)
	score(scoring, qname, bg_node_file, edge_file, bg_output_file, parameters)
	scoring = "netzcore"
	parameters={"n_iteration":5, "n_sample":100, "sampling_prefix":scoring_folder+"sampled_graph."}
	score(scoring, qname, node_file, edge_file, output_file, parameters)
	score(scoring, qname, bg_node_file, edge_file, bg_output_file, parameters)
	scoring = "netshort"
	parameters={"nd_edge_file":nd_edge_file}
	score(scoring, qname, node_file, edge_file, output_file, parameters)
	score(scoring, qname, bg_node_file, edge_file, bg_output_file, parameters)
	score_combined([output_file+".netscore", output_file+".netzcore", output_file+".netshort"], output_file+".netcombo")
	score_combined([bg_output_file+".netscore", bg_output_file+".netzcore", bg_output_file+".netshort"], bg_output_file+".netcombo")
    else:
	score(scoring_type, qname, node_file, edge_file, output_file, parameters)
	score(scoring_type, qname, bg_node_file, edge_file, bg_output_file, parameters)
    return

def score_combined(scores_file_list, output_scores_file, combination_type="standard", reverse_ranking=False):
    """
	Calculates a combined score based on normalized scores of each scoring method
    """
    node_to_scores = {}
    inf = float("Inf")
    for scores_file in scores_file_list:
	node_to_score_inner = {}
	for line in open(scores_file):
	    node, score = line.strip().split() 
	    score = float(score)
	    if inf == score:
		score = 999999 # hard coded score to correspond infinity in func. flow
	    node_to_score_inner[node] = score
	if combination_type == "standard":
	    mean, sigma = stat_utilities.calc_mean_and_sigma(node_to_score_inner.values())
	    for node, score in node_to_score_inner.iteritems():
		if sigma == 0:
		    if score-mean == 0:
			node_to_scores.setdefault(node, []).append(0)
		    else:
			node_to_scores.setdefault(node, []).append(float("inf"))
		else:
		    node_to_scores.setdefault(node, []).append((score-mean)/sigma)
	else:
	    for node, score in node_to_score_inner.iteritems():
		node_to_scores.setdefault(node, []).append(score)
    values = []
    for node, scores in node_to_scores.iteritems():
	if combination_type == "standard":
	    score = sum(scores) / len(scores)
	elif combination_type == "max":
	    score = max(scores)
	elif combination_type == "min":
	    score = min(scores)
	else:
	    raise ValueError("Unknown combination type " + combination_type)
	values.append((score, node))
    values.sort()
    min_v, max_v = min(values)[0], max(values)[0]
    f = open(output_scores_file, 'w')
    for score, node in values:
	score = (score-min_v) / (max_v-min_v)
	if reverse_ranking:
	    score = 1 - score
	f.write("%s\t%f\n" % (node, score))
    f.close()
    return

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    import numpy as np
    pvalues = np.array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = np.zeros(n)
    if correction_type == "Bonferroni":
	new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
	values.sort()
	for rank, vals in enumerate(values):
	    pvalue, i = vals
	    new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
	values.sort()
	for rank, vals in enumerate(values):
	    pvalue, i = vals
	    new_pvalues[i] = (n/(rank+1)) * pvalue
    return new_pvalues

def get_significance_among_node_scores(node_to_score, with_replacement=True, background_to_score=None, seeds=None, random_seed=None):
    from random import randint, shuffle, seed
    node_to_significance = {}
    if with_replacement:
	# 10000 times selects a node from network and checks how many of these cases 
	# the selected node has a score greater or equal to node in concern
	scores = background_to_score.values()
	size = len(scores)-1
	for node, score in node_to_score.iteritems():
	    n = 0
	    for i in xrange(10000): 
		selected = scores[randint(0,size)]
		if selected >= score:
		    n += 1
	    node_to_significance[node] = n/10000.0
    else:
	# Selects 1000 nodes and checks how many of them has a score 
	# greater or equal to node in concern
	from numpy import array
	seed(random_seed) # if None current system time is used
	scores = [ score for node, score in node_to_score.iteritems() if node not in seeds ]
	shuffle(scores)
	selected = array(scores[:1000])
	for node, score in node_to_score.iteritems():
	    n = (selected >= score).sum()
	    node_to_significance[node] = n/1000.0
    return node_to_significance

def get_node_to_description(node_mapping_file, network_file):
    network_nodes = get_nodes(network_file)
    print len(network_nodes)

    id_mapping = file_converter.get_id_to_mapped_id_mapping(node_mapping_file)
    node_to_desc = {}
    for node, vals in id_mapping.iteritems():
	selected_val = vals[0]
	for val in vals:
	    if val in network_nodes:
		selected_val = val
	node_to_desc[node] = selected_val 
    return node_to_desc

def get_nodes(file_name):
    nodes, dummy, dummy, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_name, store_edge_type = False)
    #nodes = set([ line.strip() for line in open(file_name) ])
    return nodes

def get_node_to_score(score_file):
    nodes, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = score_file, store_edge_type = False)
    return node_to_score

def get_top_nodes(score_file, selection_type="pvalue", background_score_file=None, seed_file=None, exclude_seeds=False):
    top_nodes = set() 
    node_to_score = get_node_to_score(score_file)
    if selection_type == "2sigma":
	from numpy import mean, std
	seeds = get_nodes(seed_file)
	values = []
	for node, score in node_to_score.iteritems():
	    if node not in seeds:
		values.append((score, node))
	    else: # include only seeds that are in the network
		if exclude_seeds == False:
		    top_nodes.add(node)
	m = mean(zip(*values)[0])
	s = std(zip(*values)[0])

	for score, node in values:
	    val = (score - m) / s
	    if val >= 2.0:
		top_nodes.add(node)
    elif selection_type == "pvalue":
	background_to_score = get_node_to_score(background_score_file)
	node_to_significance = get_significance_among_node_scores(node_to_score, with_replacement=True, background_to_score=background_to_score)
	pvalues = [ (val, node) for node, val in node_to_significance.iteritems() ]
	new_pvalues = correct_pvalues_for_multiple_testing(zip(*pvalues)[0])
	i = 0
	for node, val in node_to_significance.iteritems():
	    if new_pvalues[i] <= 0.05:
		top_nodes.add(node)
	    i += 1
    else:
	raise ValueError("Invalid selection type!")
    return top_nodes


if __name__ == "__main__":
    main()

