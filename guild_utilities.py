
from toolbox import network_utilities, file_converter, stat_utilities, selection_utilities
import os

def main():
    network_file = "interactions.sif"
    seed_file = "seeds.txt"
    scoring_folder = "./test/"
    executable_path = "/home/emre/arastirma/netzcore/src/scoreNetwork/scoreN"
    # Create input files for scoring
    prepare_scoring(network_file, seed_file, scoring_folder, non_seed_score=0.01, seed_score=1.0, edge_score=1.0, n_sample=100, delim=" ")
    # Run GUILD and create output files, the case for Netcombo
    run_scoring(scoring_folder, executable_path, scoring_type="netcombo")
    #run_scoring(scoring_folder, executable_path, scoring_type="netzcore", parameters={"n_iteration":5, "n_sample":100, "sampling_prefix":scoring_folder+"sampled_graph."}, qname=None)
    # Generate cross validation files
    node_scores_file = scoring_folder + "node_scores.sif"
    edge_scores_file = scoring_folder + "edge_scores_netshort.sif"

    # fill the code to get nodes, seed_to_score, edges and edge_to_score variables below
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    seeds = guild_utilities.get_nodes(seed_file)
    nodes = g.nodes()
    edges = g.edges()
    seed_to_score = dict([(node, 1) for node in seeds])
    edge_to_score = dict([((u,v), 1) for u,v in edges])

    guild_utilities.generate_cross_validation_node_score_files(nodes, seed_to_score, node_scores_file, xval = 3, default_score = 0.01, replicable = 123)

    guild_utilities.generate_cross_validation_edge_score_as_node_score_files(edges, seed_to_score, edge_to_score, edge_scores_file, xval = 3, default_score = 0.01, replicable = 123)

    # Run NetScore on these cross validation files
    guild_utilities.run_scoring(scoring_folder, executable_path, scoring_type="netscore", parameters={"n_iteration":2, "n_repetition":3}, qname=None, calculate_pvalue=True, xval=3)
    return

def prepare_scoring(network_file, seed_file, scoring_folder="./", non_seed_score=0.01, seed_score=1.0, edge_score=1.0, n_sample=100, delim=" ", name=None):
    """
	Creates input files required by GUILD executable.

	network_file: network in sif-like format where edge type is edge score (e.g., "A 0.5 B" or "A pp B")
	seed_file: seeds in text format where nodes and their scores are given (e.g., "A 0.1" or "A")
	scoring_folder: path to directory where the input/output files will be created
	non_seed_score: initial scores of non-seeds (0.01, by default)
	seed_score: initial scores of seeds (1.0, by default)
	edge_score: weight of edges, in case the values in network_file is not convertable to float  (1.0, by default)
	n_sample: number of randomly generated graphs for netzcore (100, by default)
	delim: delimiter that separates columns in input/output files (" ", by default)
	name: optional name defining the phenotype, the scoring files will created under this dir (in case of multiple phenotype analysis)
    """
    if not os.path.exists(scoring_folder):
	os.mkdir(scoring_folder)
    if name is not None:
	if not os.path.exists(scoring_folder+name):
	    os.mkdir(scoring_folder+name)
	name += os.sep
    else:
	name = ""

    # Read node info from network file (use network file as edge file)
    print "Creating edge score file"
    edge_score_file = scoring_folder + "edge_scores.sif" #network_file.split("/")[-1] + ".converted"
    if os.path.exists(edge_score_file):
	print "\tEdge score file exists, overwriting!"
    nodes, edges, dummy, edge_to_data = network_utilities.get_nodes_and_edges_from_sif_file(network_file, store_edge_type = True, delim = delim, data_to_float=False)
    edge_to_weight = create_edge_score_file(edge_score_file, edges, edge_to_data, edge_score, delim)

    # Create node file (ignore seeds that are not in the network and assign non-seed scores)
    print "Creating node score file"
    node_score_file = scoring_folder + name + "node_scores.sif" #seed_file.split("/")[-1] + ".converted"
    seed_score_file = scoring_folder + name + "seed_scores.sif"
    seeds, dummy, seed_to_data, dummy = network_utilities.get_nodes_and_edges_from_sif_file(seed_file, store_edge_type = False, delim = delim, data_to_float=True)
    if seed_to_data is None:
	seed_to_data = {}
	for seed in seeds: 
	    seed_to_data[seed] = seed_score
    node_to_data = create_node_score_file(node_score_file, seed_score_file, nodes, seeds, seed_to_data, non_seed_score, seed_score, delim)

    # Create background node file (selects k non-seeds randomly where k is the number of seeds)
    print "Creating background node score file"
    bg_node_file = scoring_folder + name + "node_scores_background.sif" #seed_file.split("/")[-1] + ".converted"
    bg_seed_file = scoring_folder + name + "seed_scores_background.sif" 
    create_background_score_file(bg_node_file, bg_seed_file, nodes, seeds, seed_to_data, non_seed_score, delim)

    # Create modified edge file using node scores for netshort
    print "Creating node score converted edge file (for netshort)"
    nd_edge_file = scoring_folder + name + "edge_scores_netshort.sif" #network_file.split("/")[-1] + ".converted_for_netshort"
    create_node_score_converted_edge_score_file(nd_edge_file, edges, edge_to_weight, node_to_data, delim)

    # Create random network files for netzcore
    print "Creating random networks (for netzcore)"
    sampling_prefix = scoring_folder + "../"  + "sampled_graph."
    if os.path.exists(sampling_prefix+"%s"%n_sample):
	print "\tSampled networks exists, skipping this step!"
    else:
	g = network_utilities.create_network_from_sif_file(network_file_in_sif = edge_score_file, use_edge_data = True, delim = delim)
	for i in xrange(1,n_sample+1):
	    g_sampled = network_utilities.randomize_graph(graph=g, randomization_type="preserve_topology_and_node_degree")
	    network_utilities.output_network_in_sif(g_sampled, sampling_prefix+"%s"%i)
    return


def create_edge_score_file(edge_score_file, edges, edge_to_data, edge_score=1.0, delim=" "):
    edge_to_weight = {}
    f = open(edge_score_file, 'w')
    for edge in edges:
	data = edge_to_data[edge]
	try:
	    score = float(data)
	except:
	    score = edge_score
	edge_to_weight[edge] = score
	f.write("%s%s%f%s%s\n" % (edge[0], delim, score, delim, edge[1]))
    f.close()
    return edge_to_weight


def create_node_score_file(node_score_file, seed_score_file, nodes, seeds, seed_to_data, non_seed_score=0.01, seed_score=1.0, delim=" "):
    f = open(node_score_file, 'w')
    f2 = open(seed_score_file, 'w')
    node_to_data = {}
    for node in nodes:
	if node in seeds:
	    if seed_to_data is not None:
		score = seed_to_data[node]
	    else:
		score = seed_score
	    f2.write("%s%s%f\n" % (node, delim, score))
	else:
	    score = non_seed_score
	node_to_data[node] = score
	f.write("%s%s%f\n" % (node, delim, score))
    f.close()
    f2.close()
    return node_to_data


def create_background_score_file(bg_node_file, bg_seed_file, nodes, seeds, seed_to_data, non_seed_score=0.01, delim=" "):
    from random import shuffle
    non_seeds = list(nodes - seeds)
    shuffle(non_seeds)
    random_seeds = set(non_seeds[:len(seeds)])
    #random_seeds = set() # in this case netzcore and netshort give the same scores for all the nodes
    f = open(bg_node_file, 'w')
    f2 = open(bg_seed_file, 'w')
    if seed_to_data is not None: seed_scores = seed_to_data.values()
    for node in nodes:
	if node in random_seeds:
	    if seed_to_data is not None:
		score = seed_scores.pop()
	    else:
		score = seed_score
	    f2.write("%s%s%f\n" % (node, delim, score))
	else:
	    score = non_seed_score
	f.write("%s%s%f\n" % (node, delim, score))
    f.close()
    f2.close()
    return


def create_node_score_converted_edge_score_file(nd_edge_file, edges, edge_to_weight, node_to_data, delim=" "):
    f = open(nd_edge_file, 'w')
    for u,v in edges:
	score_u = node_to_data[u]
	score_v = node_to_data[v]
	weight = edge_to_weight[(u, v)]
	f.write("%s%s%f%s%s\n" % (u, delim, weight*(score_u + score_v) / 2.0, delim, v))
    f.close()
    return


def run_scoring(scoring_folder, executable_path, scoring_type="netscore", parameters={"n_iteration":2, "n_repetition":3, "n_sample":100, "sampling_prefix":"./sampled_graph.", "nd_edge_file":"edge_scores_netshort.sif"}, qname=None, calculate_pvalue=True, xval=None, name=None):
    """
	Runs GUILD and creates output files.

	scoring_folder: path to directory where the input/output files will be created
	executable_path: path to GUILD executable
	scoring_type: type of the scoring method to apply, one of "netscore" | "netzcore" | "netshort" | "netcombo"
	parameters: dictionary containing scoring method specific arguments, possible keys are "n_iteration" | "n_repetition" | "n_sample" | "sampling_prefix" | "nd_edge_file"
	qname: name of the queue to be used in the cluster, one of None | "sbi" | "sbi-short" | "bigmem"
	calculate_pvalue: if True runs the scoring methods on the same network using random seeds and calculates a p-value in addition to the scores
	xval = if not None, number of folds in cross-validation
	name: optional name defining the phenotype, the scoring files will created under this dir (in case of multiple phenotype analysis)
    """

    def score(scoring_type, qname, node_file, edge_file, output_file, parameters):
	output_file += ".%s" % scoring_type
	if scoring_type == "netscore": 
	    score_command = executable_path + ' -s s -n "%s" -e "%s" -o "%s" -r %d -i %d' % (node_file, edge_file, output_file, parameters["n_repetition"], parameters["n_iteration"])
	elif scoring_type == "netzcore": 
	    score_command = executable_path + ' -s z -n "%s" -e "%s" -o "%s" -i %d -x %d -d "%s"' % (node_file, edge_file, output_file, parameters["n_iteration"], parameters["n_sample"], parameters["sampling_prefix"])
	elif scoring_type == "netshort": 
	    score_command = executable_path + ' -s d -n "%s" -e "%s" -o "%s"' % (node_file, parameters["nd_edge_file"], output_file)
	else:
	    raise ValueError("Invalid scoring type!")
	if qname is None:
	    os.system(score_command)
	else:
	    os.system("qsub -cwd -o out -e err -q %s -N %s -b y %s" % (qname, scoring_type, score_command))
	return

    if name is not None:
	name += os.sep
    else:
	name = ""

    edge_file = scoring_folder + "edge_scores.sif" 
    node_file = scoring_folder + name + "node_scores.sif" 
    seed_file = scoring_folder + name + "seed_scores.sif" 
    bg_node_file = scoring_folder + name + "node_scores_background.sif" 
    bg_seed_file = scoring_folder + name + "seed_scores_background.sif" 
    nd_edge_file = scoring_folder + name + "edge_scores_netshort.sif"
    sampling_prefix = scoring_folder + "../" + "sampled_graph."
    output_file = scoring_folder + name + "output_scores.sif"
    bg_output_file = scoring_folder + name + "output_scores_background.sif"
    if not os.path.exists(node_file) or not os.path.exists(edge_file):
	print "Input files not found!\nMake sure that you have run prepare_scoring first and that you provide the correct path."
	return
    # Run scoring algorithm
    if scoring_type == "netcombo":

	# NetScore, NetZcore, NetShort
	for scoring in ("netscore", "netzcore", "netshort"):
	    if scoring == "netscore":
		parameters={"n_repetition":3, "n_iteration":2}
	    elif scoring == "netzcore":
		parameters={"n_iteration":5, "n_sample":100, "sampling_prefix":sampling_prefix}
	    elif scoring == "netshort":
		parameters={"nd_edge_file":nd_edge_file}
	    else:
		raise ValueError("Scoring type not recognized! " + scoring)
	    score(scoring, qname, node_file, edge_file, output_file, parameters)
	    if xval is not None:
		for i in range(1, xval+1):
		    score(scoring, qname, node_file+".%d" % i, edge_file, output_file+".%d" % i, parameters)
	    if calculate_pvalue:
		score(scoring, qname, bg_node_file, edge_file, bg_output_file, parameters)
		if xval is not None:
		    for i in range(1, xval+1):
			score(scoring, qname, bg_node_file+".%d" % i, edge_file, bg_output_file+".%d" % i, parameters)

	# Netcombo
	score_combined([output_file+".netscore", output_file+".netzcore", output_file+".netshort"], output_file+".netcombo")
	if xval is not None:
	    for i in range(1, xval+1):
		score_combined(scoring, qname, node_file+".%d" % i, edge_file, output_file+".%d" % i, parameters)
	if calculate_pvalue:
	    score_combined([bg_output_file+".netscore", bg_output_file+".netzcore", bg_output_file+".netshort"], bg_output_file+".netcombo")
	    output_pvalue_file(output_file+".netcombo", bg_output_file+".netcombo", seed_file, bg_seed_file)
	    if xval is not None:
		for i in range(1, xval+1):
		    score_combined(scoring, qname, bg_node_file+".%d" % i, edge_file, bg_output_file+".%d" % i, parameters)
		    output_pvalue_file(output_file+".netcombo"+".%d" % i, bg_output_file+".netcombo"+".%d" % i, None, bg_seed_file)
    else:
	#if "sampling_prefix" not in parameters:
	#    parameters["sampling_prefix"] = sampling_prefix
	score(scoring_type, qname, node_file, edge_file, output_file, parameters)
	if xval is not None:
	    for i in range(1, xval+1):
		score(scoring_type, qname, node_file+".%d" % i, edge_file, output_file+".%d" % i, parameters)
	if calculate_pvalue:
	    score(scoring_type, qname, bg_node_file, edge_file, bg_output_file, parameters)
	    output_pvalue_file(output_file+"."+scoring_type, bg_output_file+"."+scoring_type, seed_file, bg_seed_file)
	    if xval is not None:
		for i in range(1, xval+1):
		    score(scoring_type, qname, bg_node_file+".%d" % i, edge_file, bg_output_file+".%d" % i, parameters)
		    output_pvalue_file(output_file + ".%d" % i + "."+scoring_type, bg_output_file+"."+scoring_type, None, bg_seed_file)
    return

def score_combined(scores_file_list, output_scores_file, combination_type="standard", reverse_ranking=False):
    """
	Calculates a combined score based on normalized scores of each scoring method.
	NetCombo is a special case of this method where the scores from NetShort, NetScore and NetZcore are combined.
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

def get_values_from_pvalue_file(pvalue_file, delim=" "):
    """
	Parses p-value file created using GUILD scores (e.g., using run_scoring, scoring_folder/output_scores.sif.scoring_type.pval).
    """
    f = open(pvalue_file)
    line = f.readline() # skip header
    node_to_values = {}
    for line in f.readlines():
	node, score, pval = line.strip().split(delim)
	node_to_values[node] = (float(score), float(pval))
    return node_to_values

def output_pvalue_file(score_file, background_file, seed_file=None, background_seed_file=None, delim=" "):
    """
	Calculates and outputs p-values using GUILD scores.
    """
    node_to_score = get_node_to_score(score_file)
    background_to_score = get_node_to_score(background_file)
    seed_to_score = None
    background_seed_to_score = None
    if seed_file is not None:
	seed_to_score = get_node_to_score(seed_file)
    if background_seed_file is not None:
	background_seed_to_score = get_node_to_score(background_seed_file)
	for seed in background_seed_to_score:
	    del background_to_score[seed]
    node_to_significance = get_significance_among_node_scores(node_to_score, background_to_score)
    #pvalues = node_to_significance.values()
    #new_pvalues = correct_pvalues_for_multiple_testing(pvalues) 
    values = [ (v, k) for k,v in node_to_significance.iteritems() ]
    values.sort()
    i = 0
    f = open(score_file + ".pval", 'w')
    f.write("Id%sScore%sP-value\n" % (delim, delim)) #Adjusted_P-value
    for val, node in values: #node_to_significance.iteritems():
	if seed_to_score is not None and node in seed_to_score:
	    f.write("%s%s%f%s%s\n" % (node, delim, node_to_score[node], delim, 0))
	else:
	    f.write("%s%s%f%s%s\n" % (node, delim, node_to_score[node], delim, str(val)))
	i += 1
    f.close()
    return

def output_edge_pvalue_file(network_file, score_file, background_file, seed_file=None, background_seed_file=None, delim=" "):
    """
    Calculates and outputs edge p-values using GUILD scores (node scores are first converted to edge scores and then edge p-values are calculated)
    """
    g = network_utilities.create_network_from_sif_file(network_file)
    node_to_score = get_node_to_score(score_file)
    background_to_score = get_node_to_score(background_file)
    seed_to_score = None
    background_seed_to_score = None
    #if seed_file is not None:
    #	seed_to_score = get_node_to_score(seed_file)
    if background_seed_file is not None:
	background_seed_to_score = get_node_to_score(background_seed_file)
    edge_to_score = {}
    background_edge_to_score = {}
    for u, v in g.edges():
	edge_to_score[(u,v)] = (node_to_score[u] + node_to_score[v]) / 2
	if u in background_seed_to_score or v in background_seed_to_score:
	    continue
	background_edge_to_score[(u,v)] = (background_to_score[u] + background_to_score[v]) / 2
    node_to_significance = get_significance_among_node_scores(edge_to_score, background_edge_to_score)
    values = [ (v, k) for k,v in node_to_significance.iteritems() ]
    values.sort()
    i = 0
    f = open(score_file + ".edge_pval", 'w')
    f.write("Id1%sId2%sScore%sP-value\n" % (delim, delim, delim)) 
    for val, edge in values:
	f.write("%s%s%s%s%f%s%s\n" % (edge[0], delim, edge[1], delim, edge_to_score[edge], delim, str(val)))
	i += 1
    f.close()
    return

def get_significance_among_node_scores(node_to_score, background_to_score, n_fold=10000, n_sample = 1000): 
    """
	Calculates p-values for given GUILD scores using background GUILD scores (i.e. scores calculated using random seeds).

	node_to_score: dictionary containing nodes and calculated GUILD scores
	background_to_score: dictionary containing nodes and calculated GUILD background scores
	n_fold: repeats the procedure this number of times to get p-values
	n_sample: this number of times selects a node from network and checks how many of these cases 
	the selected node has a score greater or equal to score cutoff (score bins)
    """
    from numpy import empty, array, arange, searchsorted, mean
    node_to_significance = {}
    score_cutoffs = arange(0, 1.01, 0.01)
    n_bin = len(score_cutoffs)
    folds = empty((n_fold, n_bin))
    values = background_to_score.values()
    for i, selected in enumerate(selection_utilities.get_subsamples(values, n_fold, n_sample)):
	selected = array(selected)
	bins = empty(n_bin)
	for j, score in enumerate(score_cutoffs):
	    n = (selected >= score).sum()
	    bins[j] = n/float(n_sample)
	folds[i,:] = bins
    # Average values
    folds = mean(folds, axis=0)
    #print "c(%s)" % ", ".join([ str(i) for i in folds])
    for node, score in node_to_score.iteritems():
	# Find the bin the score falls under
	index = searchsorted(score_cutoffs, score)
	if score_cutoffs[index] != score and score != 0:
	    index -= 1
	node_to_significance[node] = folds[index]
	#print node, score, folds[index]
    return node_to_significance

def get_node_to_description(node_mapping_file, network_file):
    """
	Parses tab separated node to description file.
    """
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
    """
	Parses nodes from a given file (e.g., seed file).
    """
    nodes, dummy, dummy, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_name, store_edge_type = False)
    #nodes = set([ line.strip() for line in open(file_name) ])
    return nodes

def get_node_to_score(score_file):
    """
	Parses scores from a scoring file created by GUILD.
    """
    nodes, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = score_file, store_edge_type = False)
    return node_to_score

def get_top_nodes(pvalue_file, selection_type="pvalue", seed_file=None, cutoff=None, exclude_seeds=False, delim=" "):
    """
	Selects top scoring nodes using GUILD p-value file.

	pvalue_file: GUILD p-value file
	selection_type: 2sigma | pvalue
	seed_file: file containing seeds
	cutoff: pvalue (e.g., 0.05) or sigma (e.g., 2) cutoff
	exclude_seeds: If True all the seeds given in seed file are excluded during selection
    """
    import TsvReader
    top_nodes = set() 
    reader = TsvReader.TsvReader(pvalue_file, delim=delim)
    #["Id", "Score", "P-value", "Adjusted_P-value"]
    columns, node_to_values = reader.read(fields_to_include = None)
    if selection_type == "sigma":
	if cutoff is None: cutoff = 2.0
	from numpy import mean, std
	seeds = get_nodes(seed_file)
	values = []
	for node, vals in node_to_values.iteritems():
	    try:
		# p-val file
		score, pval, adj_pval = vals[0]
	    except:
		# output scores file
		score = vals[0][0]
	    if node not in seeds:
		values.append((float(score), node))
	    else: # include only seeds that are in the network
		if exclude_seeds == False:
		    top_nodes.add(node)
	m = mean(zip(*values)[0])
	s = std(zip(*values)[0])

	for score, node in values:
	    val = (score - m) / s
	    if val >= cutoff:
		top_nodes.add(node)
    elif selection_type.startswith("pvalue"):
	if cutoff is None: cutoff = 0.05
	#nodes = empty(len(node_to_values),dtype="a16") # uses numpy, faster
        #pvalues = empty(len(node_to_values))
	#i = 0
	for node, values in node_to_values.iteritems():
	    score, pval = values[0] #, adj_pval = values[0]
	    #print node, score, pval, adj_pval
	    val = float(pval)
	    #if selection_type=="pvalue-adj":
	    #	val = float(adj_pval)
	    if val <= cutoff: 
	    	top_nodes.add(node)
	    #nodes[i] = node
	    #pvalues[i] = val 
	    #i += 1
	#top_nodes = nodes[pvalues<=cutoff]

    else:
	raise ValueError("Invalid selection type!")
    return top_nodes

def generate_cross_validation_node_score_files(nodes, seed_to_score, node_scores_file, xval = 5, default_score = 0.01, replicable = 123):
    """
	Creates GUILD input node score files for k-fold cross validation.
	2*k files are created where k of them are training files and the other k are test files.

	nodes: all the nodes in the network
	seed_to_score: dictionary containing seeds and their scores
	node_scores_file: name of the output file (".i" and ".i.test" will be appended to the training and test files for i in {1, .., k}) 
	xval: k in k-fold cross validation
	default_score: default non-seed score
	replicable: either None or an integer (that is going to be used to set the state of the random yielding in the same selection of randoms when repeated)
    """
    seeds = seed_to_score.keys()
    for k, training, test in selection_utilities.k_fold_cross_validation(seeds, xval, randomize = True, replicable = replicable):
	create_node_scores_file(nodes = nodes, node_to_score = seed_to_score, node_scores_file = node_scores_file+".%i"%k, ignored_nodes = test, default_score = default_score )
	create_node_scores_file(nodes = test, node_to_score = seed_to_score, node_scores_file = node_scores_file+".%i.test"%k, ignored_nodes = None, default_score = default_score )
    return

def generate_cross_validation_edge_score_as_node_score_files(edges, seed_to_score, edge_to_score, edge_scores_file, xval = 5, default_score = 0.01, replicable = 123):
    """
	Creates GUILD input node-converted edge score files (required by NetShort) for k-fold cross validation.
	2*k files are created where k of them are training files and the other k are test files.

	edges: edges of the network
	seed_to_score: dictionary containing seeds and their scores
	edge_to_score: dictionary containing edges and their scores
	edge_scores_file: name of the output file (".i" and ".i.test" will be appended to the training and test files for i in {1, .., k})
	xval: k in k-fold cross validation
	default_score: default non-seed score
	replicable: either None or an integer (that is going to be used to set the state of the random yielding in the same selection of randoms when repeated)
    """
    seeds = seed_to_score.keys()
    for k, training, test in selection_utilities.k_fold_cross_validation(seeds, xval, randomize = True, replicable = replicable):
	create_edge_scores_as_node_scores_file(edges = edges, node_to_score = seed_to_score, edge_to_score = edge_to_score, edge_scores_file = edge_scores_file+".%i"%k, ignored_nodes = test, default_score = default_score)
    return

def create_node_scores_file(nodes, node_to_score, node_scores_file, ignored_nodes = None, default_score = 0, delim=" "):
    """
	Creates node association/relevance score file for the nodes in the network.

	ignored_nodes: nodes in this set are ignored (assigned default_score)
    """
    f = open(node_scores_file, 'w')
    for v in nodes:
	if ignored_nodes is not None and v in ignored_nodes:
	    score_v = default_score
	else:
	    if node_to_score.has_key(v):
		score_v = node_to_score[v]
	    else:
		score_v = default_score
	f.write("%s%s%f\n" % (v, delim, score_v))
    f.close()
    return 

def create_edge_scores_as_node_scores_file(edges, node_to_score, edge_to_score, edge_scores_file, ignored_nodes = None, default_score = 0, delim=" "):
    """
	Creates edge score file from node association scores required by NetShort.
    """
    f = open(edge_scores_file, 'w')
    for u,v in edges:
	if ignored_nodes is not None and u in ignored_nodes:
	    score_u = default_score
	else:
	    if node_to_score.has_key(u):
		score_u = node_to_score[u]
	    else:
		score_u = default_score
	if ignored_nodes is not None and v in ignored_nodes:
	    score_v = default_score
	else:
	    if node_to_score.has_key(v):
		score_v = node_to_score[v]
	    else:
		score_v = default_score
	weight = 1  
	if (u,v) in edge_to_score:
	    weight = edge_to_score[(u,v)]
	elif (v,u) in edge_to_score:
	    weight = edge_to_score[(v,u)]
	f.write("%s%s%f%s%s\n" % (u, delim, weight*(score_u + score_v) / 2.0, delim, v))
    f.close()
    return

if __name__ == "__main__":
    main()

