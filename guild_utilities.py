
from toolbox import network_utilities, file_converter

def main():
    return

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

def get_top_nodes(score_file, seed_file):
    from numpy import mean, std

    top_nodes = set() 
    node_to_score = get_node_to_score(score_file)
    seeds = get_nodes(seed_file)
    values = []
    for node, score in node_to_score.iteritems():
	if node not in seeds:
	    values.append((score, node))
	else: # include only seeds that are in the network
	    top_nodes.add(node)
    m = mean(zip(*values)[0])
    s = std(zip(*values)[0])

    for score, node in values:
	val = (score - m) / s
	if val >= 2.0:
	    top_nodes.add(node)
    return top_nodes


if __name__ == "__main__":
    main()

