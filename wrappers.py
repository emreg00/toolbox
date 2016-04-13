#######################################################################
# Recipies / wrapper functions using toolbox methods for disease,
# drug and network analysis
# e.g. 10/2015
#######################################################################
import network_utilities, parse_msigdb, stat_utilities, dict_utilities
import parse_uniprot, parse_ncbi
import csv, numpy, os, cPickle

##### Id mapping related #####

def get_uniprot_to_geneid(uniprot_file, uniprot_ids):
    uniprot_to_geneid = parse_uniprot.get_uniprot_to_geneid(uniprot_file, uniprot_ids)
    return uniprot_to_geneid


def get_geneid_symbol_mapping(mapping_file):
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    return geneid_to_names, name_to_geneid


##### Network related #####

def get_network(network_file, only_lcc):
    network = network_utilities.create_network_from_sif_file(network_file, use_edge_data = False, delim = None, include_unconnected=True)
    #print len(network.nodes()), len(network.edges())
    if only_lcc:
	components = network_utilities.get_connected_components(network, False)
	network = network_utilities.get_subgraph(network, components[0])
	#print len(network.nodes()), len(network.edges())
    return network


##### Gene expression related #####

def get_expression_info(gexp_file, process=None, delim=',', quote='"', dump_file=None):
    """
    To get gene expression info
    process: a set(["log2", "z", "abs"]) or None
    """
    if dump_file is not None and os.path.exists(dump_file):
	gexp, gene_to_idx, cell_line_to_idx = cPickle.load(open(dump_file))
	return gexp, gene_to_idx, cell_line_to_idx
    #gene_to_values = {}
    f = open(gexp_file)
    reader = csv.reader(f, delimiter=delim, quotechar=quote)
    header = reader.next()
    #print len(header), header
    header = header[1:]
    cell_line_to_idx = dict([ (cell_line, i) for i, cell_line in enumerate(header) ])
    gene_to_idx = {}
    values_arr = []
    for i, row in enumerate(reader):
	gene = row[0]
	values = map(float, row[1:])
	#gene_to_values[gene] = values
	gene_to_idx[gene] = i
	values_arr.append(values)
    f.close()
    gexp = numpy.array(values_arr)
    if process is not None:
	if "log2" in process:
	    gexp = numpy.log2(gexp)
	if "z" in process:
	    gexp = (gexp - gexp.mean(axis=1)[:, numpy.newaxis]) / gexp.std(axis=1, ddof=1)[:, numpy.newaxis]
	if "abs" in process:
	    gexp = numpy.abs(gexp)
	#print gexp.shape, gexp_norm.shape
	#print gexp[0,0], gexp_norm[0,0]
	#return gene_to_values, cell_line_to_idx
    if dump_file is not None:
	values = gexp, gene_to_idx, cell_line_to_idx
	cPickle.dump(values, open(dump_file, 'w')) 
    return gexp, gene_to_idx, cell_line_to_idx


##### Pathway related #####

def get_pathway_info(pathway_file, prefix=None, nodes=None):
    """
    nodes to filter geneids that are not in the network
    prefix: kegg | reactome | biocarta
    """
    pathway_to_geneids, geneid_to_pathways = parse_msigdb.get_msigdb_info(pathway_file, prefix)
    if nodes is not None:
	pathway_to_geneids_mod = {}
	for pathway, geneids in pathway_to_geneids.iteritems():
	    geneids_mod = geneids & nodes
	    if len(geneids_mod) > 0:
		pathway_to_geneids_mod[pathway] = geneids_mod 
	pathway_to_geneids = pathway_to_geneids_mod
    return pathway_to_geneids


def get_diseasome_genes(diseasome_file, nodes=None):
    disease_to_genes = {}
    disease_to_category = {}
    for line in open(diseasome_file):
	words = line.strip("\n").split("\t")
	disease = words[1].strip('"')
	category = words[0]
	genes = set(words[2:])
	if nodes is not None:
	    genes &= nodes
	    if len(genes) == 0:
		continue
	disease_to_genes[disease] = genes
	disease_to_category[disease] = category
    return disease_to_genes, disease_to_category


##### Statistics related #####

def overlap_significance(geneids1, geneids2, nodes):
    n1, n2 = len(geneids1), len(geneids2)
    n = len(geneids1 & geneids2)
    N = len(nodes)
    pval = stat_utilities.hypergeometric_test_numeric(n, n1, N, n2)
    return n, n1, n2, pval


##### Proximity related #####
def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, n_random=1000, min_bin_size=100, seed=452456):
    #distance = "closest"
    #lengths = network_utilities.get_shortest_path_lengths(network, "../data/toy.sif.pcl")
    #d = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
    nodes_network = set(network.nodes())
    if len(set(nodes_from) & nodes_network) == 0 or len(set(nodes_to) & nodes_network) == 0:
	return None # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to)
    if nodes_from_random is None:
	nodes_from_random = get_random_nodes(nodes_from, network, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
	nodes_to_random = get_random_nodes(nodes_to, network, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(n_random)
    for i, values_random in enumerate(random_values_list):
	nodes_from, nodes_to = values_random
	#values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
	values[i] = calculate_closest_distance(network, nodes_from, nodes_to)
    #pval = float(sum(values <= d)) / len(values)
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
	z = 0.0
    else:
	z = (d - m) / s
    return d, z, (m, s) #(z, pval)


def calculate_closest_distance(network, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from:
	values = []
	for node_to in nodes_to:
	    val = network_utilities.get_shortest_path_length_between(network, node_from, node_to)
	    values.append(val)
	d = min(values)
	#print d,
	values_outer.append(d)
    d = numpy.mean(values_outer)
    #print d
    return d


def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
	# Get degree bins of the network
	bins = network_utilities.get_degree_binning(network, min_bin_size) 
    nodes_random = network_utilities.pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed) 
    return nodes_random


