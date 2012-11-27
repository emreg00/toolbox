import network_utilities

def main():
    module_file = "modules.txt"
    network_file = "../data/interactions.sif"
    g = create_network_from_sif_file(network_file)
    modules = get_modules_of_graph(g, "mcl", output_file=module_file, inflation=1.7)
    # modules = get_modules_from_file(module_file)
    return

def create_network_from_sif_file(network_file, **kwargs):
    return network_utilities.create_network_from_sif_file(network_file, **kwargs)

#! Below two function need to be checked
def get_seeds_from_node_scores_file(node_scores_file, default_non_seed_score):
    nodes, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)
    seeds = set()
    for node in initial_node_to_score:
	if initial_node_to_score[node] > default_non_seed_score:
	    seeds.add(node)
    return seeds, nodes


def score_mcl(node_scores_file, network_file, output_scores_file, module_file, default_non_seed_score):
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    #modules = get_modules_of_graph(g, "mcl", inflation=2) # if edge weight based clustering is desired
    seeds, nodes = get_seeds_from_node_scores_file(node_scores_file, default_non_seed_score)
    modules = get_modules_from_file(module_file)
    f = open(output_scores_file, 'w')
    node_to_score = {}
    #selected = set()
    for module in modules:
	module = set(module)
	#common = module&seeds
	#if 100*float(len(common))/len(module) > threshold:
	    #selected |= module
	#score = float(len(common))/len(module)
	#n = len(module)-len(common)
	#if n == 0:
	#    continue
	#score = 1.0/n
	for node in module:
	    #node_to_score[node] = score
	    neighbors = set(g.neighbors(node))
	    common = neighbors & module
	    if node in common:
		common.remove(node)
	    #if len(common) == 0:
	    #	continue
	    score = float(len(common&seeds)) / len(module)
	    node_to_score[node] = score
    for node in nodes:
	if node in node_to_score:
	    f.write("%s\t%f\n" % (node, node_to_score[node]))
	else:
	    f.write("%s\t0.0\n" % node)
    f.close()
    return

def get_modules_from_file(output_file):
    f = open(output_file)
    modules = []
    for line in f:
	words = line.strip().split("\t")
	modules.append(words)
    f.close()
    return modules

def get_modules_of_graph(sub_graph, module_detection_type, output_file, inflation=1.7):
    if module_detection_type == "connected":
	import network_utilities
	modules = network_utilities.get_connected_components(sub_graph, return_as_graph_list=True)
    elif module_detection_type == "mcl":
	from os import system
	f = open(output_file + ".mcl", 'w')
	nodes = set()
	for node1, node2, data in sub_graph.edges(data=True):
	    nodes.add(node1)
	    nodes.add(node2)
	    if 'w' in data:
		data = str(data['w'])
	    else:
		data = "-"
	    f.write("%s\t%s\t%s\n" % (node1, node2, data))
	for node in sub_graph.nodes():
	    if node not in nodes:
		f.write("%s\n" % node)
	f.close()
	# Optimum inflation parameter was 1.7-1.8 in a recent comparison paper
	system("mcl %s --abc -I %f -o %s 2>> %s" % (output_file + ".mcl", inflation, output_file, output_file + ".err"))
	modules = get_modules_from_file(output_file)
    else:
	raise ValueError("Unrecognized module detection type")
    #print len(modules), map(len, modules)
    return modules

if __name__ == "__main__":
    main()


