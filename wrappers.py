#######################################################################
# Recipies / wrapper functions using toolbox methods for disease,
# drug and network analysis
# e.g. 10/2015
#######################################################################
import network_utilities, parse_msigdb, stat_utilities, dict_utilities, TsvReader
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
	network_lcc_file = network_file + ".lcc"
	if not os.path.exists(network_lcc_file ):
	    f = open(network_lcc_file, 'w')
	    for u,v in network.edges():
		f.write("%s 1 %s\n" % (u, v))
	    f.close()
    return network


##### Gene expression related #####

def get_expression_info(gexp_file, process=None, delim=',', quote='"', R_header=False, dump_file=None):
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
    if R_header == False:
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


def get_de_genes(file_name, cutoff_adj = 0.05, cutoff_fc=1.5, n_top=None, id_type = "GeneID"):
    """
    For parsing DE file generated using R PEPPER package
    """
    fields_to_include = [id_type, "P.Value", "logFC", "adj.P.Val"]
    parser = TsvReader.TsvReader(file_name, delim="\t", inner_delim=None)
    header_to_idx, id_to_values = parser.read(fields_to_include, keys_to_include=None, merge_inner_values=False)
    if "" in id_to_values:
	del id_to_values[""]
    #print len(id_to_values)
    #gene = "10458"
    #if gene in id_to_values:
    #    print id_to_values[gene] 
    genes = set()
    genes_all = set()
    values_gene = []
    for gene, values in id_to_values.iteritems():
	include = False
	for val in values:
	    pval = val[header_to_idx["adj.p.val"]] # "p.value"]]
	    if pval == "NA":
		continue
	    fc = float(val[header_to_idx["logfc"]])
	    if float(pval) <= cutoff_adj:
		if abs(fc) >= cutoff_fc: # fc >= cutoff_fc: 
		    include = True
		if n_top is not None:
		    values_gene.append((abs(fc), gene))
	for word in gene.split("///"):
	    word = word.strip()
	    if word == "---":
		continue
	    if include:
		genes.add(word)
	    else:
		genes_all.add(word)
    if n_top is not None:
	values_gene.sort()
	genes = set([ word.strip() for fc, gene in values_gene[-n_top:] for word in gene.split("///") ])
    return genes, genes_all


def get_z_genes(file_name, cutoff_z = 2):
    """
    For parsing DE-Z file generated using R PEPPER package
    """
    gexp, gene_to_idx, cell_line_to_idx = get_expression_info(file_name, process=None, delim='\t', R_header=True) #, quote='"', dump_file=None)
    genes = gene_to_idx.items()
    genes.sort(key=lambda x: x[1])
    genes = numpy.array(zip(*genes)[0])
    sample_to_genes = {}
    for cell_line, idx in cell_line_to_idx.iteritems():
	indices = numpy.abs(gexp[:,idx]) > cutoff_z 
	sample_to_genes[cell_line] = genes[indices]
	#if cell_line in ["GSM734834", "GSM734833"]: 
	#    print cell_line, len(genes[indices]), genes[indices] 
    return sample_to_genes


def get_sample_mapping(file_name, labels_case, labels_control=None):
    f = open(file_name)
    labels_case = set(labels_case)
    if labels_control is not None:
	labels_control = set(labels_control)
    samples_case = []
    samples_control = []
    for line in f:
	sample, label = line.strip("\n").split("\t")
	label = label.strip()
	if label in labels_case:
	    samples_case.append(sample)
	else:
	    if labels_control is None or label in labels_control:
		samples_control.append(sample)
    return samples_case, samples_control


##### Disease, pathway, comorbidity, symptom info related #####

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


def get_comorbidity_info(correlation_type="phi", only_significant=False):
    """
    Parse HuDiNe data
    correlation_type: phi (pearson correlation) | RR (favors rare disease pairs)
    """
    comorbidity_file = CONFIG.get("comorbidity_file")
    f = open(comorbidity_file)
    header_to_idx = dict((word, i) for i, word in enumerate(f.readline().strip().split("\t")))
    disease_to_disease_comorbidity = {}
    for line in f:
	words = line.strip().split("\t")
	disease1, disease2 = words[:2]
	significance = words[header_to_idx["sign_"+correlation_type]]
	if only_significant and significance == "0":
	    continue
	val = float(words[header_to_idx[correlation_type]])
	disease_to_disease_comorbidity.setdefault(disease1, {})[disease2] = (val, significance)
	disease_to_disease_comorbidity.setdefault(disease2, {})[disease1] = (val, significance)
    f.close()
    return disease_to_disease_comorbidity


def get_symptom_info():
    disease_to_symptoms = {}
    symptom_to_diseases = {}
    f = open(CONFIG.get("symptom_file"))
    for line in f:
	words = line.strip("\n").split("\t")
	symptom, disease, n, score = words
	symptom = symptom.lower()
	disease = disease.lower()
	#if float(score) < 3:
	#    continue
	disease_to_symptoms.setdefault(disease, set()).add(symptom)
	symptom_to_diseases.setdefault(symptom, set()).add(disease)
    return disease_to_symptoms, symptom_to_diseases


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


### GUILD related ###

def create_node_file(node_to_score, nodes, node_file, background_score = 0.01):
    """
    Simplified method for creating guild node score files
    """
    f = open(node_file, 'w')
    for node in nodes:
	if node in node_to_score:
	    score = node_to_score[node]
	else:
	    score = background_score
	f.write("%s %f\n" % (node, score))
    f.close()
    return


def run_guild(phenotype, node_to_score, network_nodes, network_file, output_dir, executable_path = None, background_score = 0.01, qname=None): 
    # Create node file
    node_file = "%s%s.node" % (output_dir, phenotype)
    create_node_file(node_to_score, network_nodes, node_file, background_score)
    output_file = "%s%s.ns" % (output_dir, phenotype)
    n_repetition = 3 
    n_iteration = 2 
    # Get and run the GUILD command
    #print strftime("%H:%M:%S - %d %b %Y") #, score_command
    score_command = ' -s s -n "%s" -e "%s" -o "%s" -r %d -i %d' % (node_file, network_file, output_file, n_repetition, n_iteration)
    if qname is None:
	if executable_path is None:
	    executable_path = "guild" # assuming accessible guild executable
	score_command = executable_path + score_command
	os.system(score_command)
    else:
	#os.system("qsub -cwd -o out -e err -q %s -N %s -b y %s" % (qname, scoring_type, score_command))
	#print "qsub -cwd -o out -e err -q %s -N guild_%s -b y %s" % (qname, drug, score_command)
	print "%s" % (score_command.replace('"', ''))
    return


### Functional enrichment related ###

def check_functional_enrichment(id_list, background_id_list = None, id_type = "genesymbol", evidences = None, out_file_name = None):
    """
    evidences = ['EXP', 'IDA', 'IEP', 'IGI', 'IMP', 'ISA', 'ISM', 'ISO', 'ISS']
    evidences = None corresponds to ['EXP', 'IC', 'IDA', 'IEA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NAS', 'RCA', 'TAS']
    for custom associations: association = [["GO:0006509", "351"], ["GO:0048167", "348", "5663", "5664", "23621"], ["GO:0097458", "1005", "1006", "1007"], ["GO:0048487", "1", "2", "351"], ["GO:0048488", map(str, range(1000,2000))]]
    """
    if out_file_name is not None:
	f_output = open(out_file_name, 'w').write
    else:
	from stdout import sys
	f_output = stdout.write
    return functional_enrichment.check_functional_enrichment(id_list, background_id_list, id_type, f_output, tex_format = False, support = evidences, associations = None) 



