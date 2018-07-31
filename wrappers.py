#######################################################################
# Recipies / wrapper functions using toolbox methods for disease,
# drug and network analysis
# e.g. 10/2015
#######################################################################
import network_utilities, stat_utilities, dict_utilities, text_utilities
import TsvReader, functional_enrichment
import parse_umls, parse_msigdb 
import parse_uniprot, parse_ncbi
import parse_do, parse_medic, parse_disgenet
import parse_drugbank, parse_medi, parse_hetionet
import csv, numpy, os, cPickle
import random
try:
    from toolbox.external.diamond import diamond
except:
    print "DIAMOnD not found and thus will not be available!"

##### Id mapping related #####

def get_mapping(file_name, from_column = None, to_column = None, delim=None, one_to_one=True):
    """
    Assumes header, maps first column to the second
    """
    key_to_value = TsvReader.get_from_to_mapping(file_name, from_column = from_column, to_column = to_column, delim=delim, inner_delim = None, filter_column = None, exclude_value = None, include_value = None, one_to_one=one_to_one) 
    return key_to_value

def convert_to_geneid(file_name, id_type, id_mapping_file):
    """
    Expects a file where each line is a gene name / uniprot id
    id_type: symbol | uniprot
    """
    genes = [ line.strip("\n") for line in open(file_name) ]
    if id_type == "symbol":
	geneid_to_name, name_to_geneid = get_geneid_symbol_mapping(id_mapping_file)
    elif id_type == "uniprot":
	name_to_geneid = get_uniprot_to_id(id_mapping_file, uniprot_ids=genes, only_min=True, key_function=int)
    else:
	raise ValueError("Uknown id type: %s" % id_type)
    geneids = set([ name_to_geneid[gene] for gene in genes if gene in name_to_geneid ])
    genes_non = set([ gene for gene in genes if gene not in name_to_geneid ])
    print "Not found genes:", genes_non
    return geneids


def get_uniprot_to_id(uniprot_file, uniprot_ids=None, only_min=True, key_function=int):
    """
    uniprot_file = %(data_dir)s/uniprot/idmapping.tab or idmapping.tab.symbol or idmapping.tab.mouse
    Can be used to convert to geneids (key_function=int) as well as symbols (key_function=len) depending on the input file
    """
    uniprot_to_gene = parse_uniprot.get_uniprot_to_geneid(uniprot_file, uniprot_ids, only_min, key_function)
    return uniprot_to_gene


def get_geneid_symbol_mapping(mapping_file):
    """
    id_mapping_file = %(data_dir)s/ncbi/geneid_to_symbol.txt
    """
    geneid_to_name, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    return geneid_to_name, name_to_geneid


def get_mesh_id_mapping(desc_file, rel_file, dump_file = None):
    """ 
    Get all concept id - mesh id mapping (also gets entry names in addition to main header)
    """
    #dump_file = CONFIG.get("umls_dir") + "/mapping.pcl"
    mesh_id_to_name, concept_id_to_mesh_id, mesh_id_to_name_with_synonyms = parse_umls.get_mesh_id_mapping(desc_file, rel_file, dump_file = dump_file)
    return mesh_id_to_name, concept_id_to_mesh_id, mesh_id_to_name_with_synonyms


def get_mesh_disease_ontology(desc_file, rel_file, dump_file = None):
    #dump_file = CONFIG.get("umls_dir") + "/ontology.pcl"
    g = parse_umls.get_mesh_disease_ontology(desc_file, rel_file, dump_file = dump_file)
    return g


def get_medic_mesh_id_mapping(medic_file):
    name_to_id, id_to_mesh_ids = get_medic_mesh_id_mapping(medic_file)
    return name_to_id, id_to_mesh_ids 


def get_mesh_disease_category_mapping(desc_file, rel_file, dump_file = None):
    #dump_file = CONFIG.get("umls_dir") + "/ontology.pcl"
    mesh_id_to_top_ids = parse_umls.get_mesh_id_to_disease_category(desc_file, rel_file, dump_file)
    mesh_id_to_name, concept_id_to_mesh_id, mesh_id_to_name_with_synonyms = get_mesh_id_mapping(desc_file, rel_file, dump_file)
    mesh_name_to_parents = {}
    for child, parents in mesh_id_to_top_ids.iteritems():
	name_child = mesh_id_to_name[concept_id_to_mesh_id[child]]
	values = []
	for parent in parents:
	    name_parent = mesh_id_to_name[concept_id_to_mesh_id[parent]]
	    values.append(name_parent.lower())
	values.sort()
	mesh_name_to_parents[name_child.lower()] = values
    #print mesh_name_to_parents["asthma"], mesh_name_to_parents["psoriasis"]
    return mesh_name_to_parents


def get_icd_to_mesh_ids(disease_ontology_file, id_type="ICD9CM"):
    icd_to_mesh_ids = parse_do.get_icd_to_mesh_ids(disease_ontology_file, id_type)
    return icd_to_mesh_ids


def get_do_to_mesh_ids(disease_ontology_file):
    name_to_do_id, do_to_mesh_ids, mesh_id_to_type_to_ids = parse_do.get_do_mesh_id_mapping(disease_ontology_file)
    return do_to_mesh_ids


def get_homology_mapping(homologene_file, tax_id="10090", from_tax_id="9606", symbol_type="geneid"):
    """
    symbol_type: geneid | symbol
    """
    geneid_to_geneid, group_to_taxid_to_geneid = parse_ncbi.get_homology_mapping(homologene_file, tax_id, from_tax_id=from_tax_id, symbol_type=symbol_type)
    return geneid_to_geneid 


##### Network related #####

def get_network(network_file, only_lcc):
    network = network_utilities.create_network_from_sif_file(network_file, use_edge_data = False, delim = None, include_unconnected=True)
    #print len(network.nodes()), len(network.edges())
    if only_lcc and not network_file.endswith(".lcc"):
	print "Shrinking network to its LCC", len(network.nodes()), len(network.edges())
	components = network_utilities.get_connected_components(network, False)
	network = network_utilities.get_subgraph(network, components[0])
	print "Final shape:", len(network.nodes()), len(network.edges())
	#print len(network.nodes()), len(network.edges())
	network_lcc_file = network_file + ".lcc"
	if not os.path.exists(network_lcc_file ):
	    f = open(network_lcc_file, 'w')
	    for u,v in network.edges():
		f.write("%s 1 %s\n" % (u, v))
	    f.close()
    return network


def create_functional_network(links_file, mapping_file, cutoff = 900): 
    #string_dir = CONFIG.get("string_dir") + "/"
    #links_file = string_dir + CONFIG.get("string_links_file")
    #mapping_file = string_dir + CONFIG.get("string_mapping_file")
    output_file = CONFIG.get("network_file")
    parse_string.get_interactions(links_file, mapping_file, output_file, cutoff) #, include_score=True) 
    #network = get_network()
    #print len(network.nodes()), len(network.edges())
    return


def calculate_lcc_significance(network, nodes, nodes_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456):
    # Degree matching problematic for small bin sizes
    #if bins is None and nodes_random is None:
    #	bins = network_utilities.get_degree_binning(network, min_bin_size) 
    random.seed(seed)
    if nodes_random is None:
	network_nodes = list(network.nodes())
	#nodes_random = get_random_nodes(nodes, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
	nodes_random = []
	for i in xrange(n_random):
	    random.shuffle(network_nodes)
	    nodes_random.append(network_nodes[:len(nodes)])
    network_sub = network.subgraph(nodes)
    component_nodes = network_utilities.get_connected_components(network_sub, False)
    #print component_nodes 
    d = len(component_nodes[0])
    values = numpy.empty(len(nodes_random)) 
    for i, nodes in enumerate(nodes_random):
	network_sub = network.subgraph(nodes)
	component_nodes = network_utilities.get_connected_components(network_sub, False)[0]
	values[i] = len(component_nodes)
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
	z = 0.0
    else:
	z = (d - m) / s
    return d, z, (m, s), values


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
	#if "na.rm" in process:
	#    idx = numpy.where(numpy.isnan(a)) # need to remove rows with NAs
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
    genes_up, genes_down = set(), set()
    values_gene = []
    for gene, values in id_to_values.iteritems():
	include = False
	positive = False
	for val in values:
	    pval = val[header_to_idx["adj.p.val"]] # "p.value"]]
	    if pval == "NA":
		continue
	    fc = float(val[header_to_idx["logfc"]])
	    if float(pval) <= cutoff_adj:
		if abs(fc) >= cutoff_fc: # fc >= cutoff_fc: 
		    include = True
		if fc >= 0: 
		    positive = True
		if n_top is not None:
		    values_gene.append((abs(fc), gene))
	for word in gene.split("///"):
	    word = word.strip()
	    if word == "---":
		continue
	    if include:
		genes.add(word)
		if positive:
		    genes_up.add(word)
		else:
		    genes_down.add(word)
	    else:
		genes_all.add(word)
    if n_top is not None:
	values_gene.sort()
	genes = set([ word.strip() for fc, gene in values_gene[-n_top:] for word in gene.split("///") ])
    return genes, genes_all, genes_up, genes_down


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

def get_pathway_info(pathway_file, prefix=None, nodes=None, max_pathway_size=None, inner_delim=None):
    """
    Assumes a tab separated file containing pathway name, link, geneids
    nodes to filter geneids that are not in the network
    prefix: kegg | reactome | biocarta
    inner_delim: None for tab separated geneids, " " for space separated geneids
    """
    pathway_to_geneids, geneid_to_pathways = parse_msigdb.get_msigdb_info(pathway_file, prefix, inner_delim=inner_delim)
    if nodes is not None or max_pathway_size is not None:
	pathway_to_geneids_mod = {}
	for pathway, geneids in pathway_to_geneids.iteritems():
	    if max_pathway_size is not None:
		if len(geneids) > max_pathway_size:
		    continue
	    if nodes is not None:
		geneids &= nodes
		if len(geneids) == 0:
		    continue
	    pathway_to_geneids_mod[pathway] = geneids
	pathway_to_geneids = pathway_to_geneids_mod
    return pathway_to_geneids


def get_diseasome_genes(diseasome_file, nodes=None, network=None):
    """
    If nodes is not None, keep only nodes in the network
    If network is not None, keep only LCC
    """
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
	if network is not None:
	    network_sub = network.subgraph(genes)
	    genes = network_utilities.get_connected_components(network_sub, False)[0]
	disease_to_genes[disease] = genes
	disease_to_category[disease] = category
    return disease_to_genes, disease_to_category


def get_disgenet_genes(file_name):
    disease_to_genes, disease_to_sources, cui_to_disease = parse_disgenet.get_disgenet_genes(file_name)
    return disease_to_genes, disease_to_sources, cui_to_disease 


def get_comorbidity_info(comorbidity_file, disease_ontology_file, mesh_dump, correlation_type="RR", only_significant=False):
    """
    Parse HuDiNe data from AllNet3 and map ICD9 to MeSH using DO
    correlation_type: phi (pearson correlation) | RR (favors rare disease pairs)
    """
    icd_to_mesh_ids = get_icd_to_mesh_ids(disease_ontology_file, id_type="ICD9CM")
    #print len(icd_to_mesh_ids), icd_to_mesh_ids.items()[:5]
    #print [ icd_to_mesh_ids[val] for val in ("289", "502", "578", "579", "542", "543")]
    mesh_id_to_name, concept_id_to_mesh_id, mesh_id_to_name_with_synonyms = get_mesh_id_mapping(None, None, dump_file = mesh_dump)
    #print mesh_id_to_name.items()[:5]
    f = open(comorbidity_file)
    header_to_idx = dict((word, i) for i, word in enumerate(f.readline().strip().split("\t")))
    disease_to_disease_comorbidity = {}
    for line in f:
	words = line.strip().split("\t")
	icd1, icd2 = words[:2]
	if only_significant and significance == "0":
	    continue
	if icd1 not in icd_to_mesh_ids or icd2 not in icd_to_mesh_ids:
	    #print "Not in DO mapping:", icd1, icd2
	    continue
	val = float(words[header_to_idx[correlation_type]]) # idx:5
	for mesh1 in icd_to_mesh_ids[icd1]:
	    if mesh1 not in mesh_id_to_name:
		#print "Not in name mapping:", mesh1
		continue
	    disease1 = mesh_id_to_name[mesh1].lower()
	    for mesh2 in icd_to_mesh_ids[icd2]:
		if mesh2 not in mesh_id_to_name:
		    #print "Not in name mapping:", mesh2
		    continue
		disease2 = mesh_id_to_name[mesh2].lower()
		disease1_mod, disease2_mod = sorted((disease1, disease2))
		d = disease_to_disease_comorbidity.setdefault(disease1, {})
		if disease2 in d:
		    if d[disease2] > val: # skip if the existing comorbidity value is higher
			continue
		d[disease2] = val 
		d = disease_to_disease_comorbidity.setdefault(disease2, {})
		d[disease1] = val
		#print icd1, mesh1, disease1, icd2, mesh2, disease2, val 
    #print len(disease_to_disease_comorbidity), disease_to_disease_comorbidity.values()[0].items()[:5]
    return disease_to_disease_comorbidity 
    # Parse HuDiNe data from potentially buggy comorbidity_new.tsv
    #comorbidity_file = CONFIG.get("comorbidity_file")
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


def get_symptom_info(symptom_file, tfidf_cutoff=None):
    """
    Parse Zhou et al supplementary s4. A cutoff of 3.5 is likely to filter spurious associations.
    """
    disease_to_symptoms = {}
    symptom_to_diseases = {}
    disease_to_symptom_to_score = {}
    #symptom_file = CONFIG.get("symptom_file")
    f = open(symptom_file)
    f.readline()
    for line in f:
	words = line.strip("\n").split("\t")
	symptom, disease, n, score = words
	symptom = symptom.lower()
	disease = disease.lower()
	if tfidf_cutoff is not None and not float(score) >= tfidf_cutoff:
	    continue
	disease_to_symptoms.setdefault(disease, set()).add(symptom)
	symptom_to_diseases.setdefault(symptom, set()).add(disease)
	d = disease_to_symptom_to_score.setdefault(disease, {})
	d[symptom] = float(score)
    return disease_to_symptoms, symptom_to_diseases, disease_to_symptom_to_score


##### Drug related info #####

def get_drugbank(drugbank_file):
    dump_file = drugbank_file + ".pcl"
    if os.path.exists(dump_file):
	parser = cPickle.load(open(dump_file))
    else:
	parser = parse_drugbank.DrugBankXMLParser(drugbank_file)
	parser.parse()
	cPickle.dump(parser, open(dump_file, 'w'))
    return parser


def get_medi_indications(medi_file, drugbank_file, mesh_dump, disease_ontology_file, only_hps=True):
    dump_file = medi_file + ".pcl"
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    parser = get_drugbank(drugbank_file)
    name_to_drug, synonym_to_drug = parser.get_synonyms(selected_drugs=None, only_synonyms=False)
    name_to_icd_and_confidences = parse_medi.get_medi_mapping(medi_file)
    mesh_id_to_name, concept_id_to_mesh_id, mesh_id_to_name_with_synonyms = get_mesh_id_mapping(None, None, dump_file = mesh_dump)
    icd_to_mesh_ids = get_icd_to_mesh_ids(disease_ontology_file, id_type="ICD9CM")
    drug_to_indications = parse_medi.get_drug_disease_mapping(name_to_icd_and_confidences, name_to_drug, synonym_to_drug, icd_to_mesh_ids, mesh_id_to_name, dump_file = None) 
    drug_to_diseases = {}
    for drug, values in drug_to_indications.iteritems():
	for phenotype, dui, val in values:
    	    if only_hps and val <= 0.5:
    		continue
	    drug_to_diseases.setdefault(drug, set()).add(phenotype)
    # Drug name to disease name textual mapping
    #name_to_indication_and_confidences = parse_medi.get_medi_mapping(medi_file, textual_indication=True)
    #drug_to_diseases = {}
    #for name, values in name_to_indication_and_confidences.iteritems():
    #	for indication, confidence in values:
    #	    if only_hps and confidence <= 0.5:
    #		continue
    #	    drug_to_diseases.setdefault(name, set()).add(indication.lower())
    # Get disease to name mapping
    #phenotype_to_mesh_id = dict((name, mesh_id) for mesh_id, name in mesh_id_to_name.iteritems())
    #disease_to_drugs = parse_medi.get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id)
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases


def get_hetionet_indications(hetionet_file, mesh_dump, disease_ontology_file):
    dump_file = hetionet_file + ".pcl"
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_do_ids = parse_hetionet.get_hetionet_mapping(hetionet_file, metaedge="CtD")
    do_to_mesh_ids = get_do_to_mesh_ids(disease_ontology_file)
    mesh_id_to_name, concept_id_to_mesh_id, mesh_id_to_name_with_synonyms = get_mesh_id_mapping(None, None, dump_file = mesh_dump)
    drug_to_indications = parse_hetionet.get_drug_disease_mapping(drug_to_do_ids, do_to_mesh_ids, mesh_id_to_name, dump_file = None)
    drug_to_diseases = {}
    for drug, values in drug_to_indications.iteritems():
	for phenotype, dui, val in values:
	    drug_to_diseases.setdefault(drug, set()).add(phenotype)
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases 


##### Statistics related #####

def overlap_significance(geneids1, geneids2, nodes, method="hyper"):
    """
    method: hyper(geometric) | fishers (two-sided version of hypergeometric) | jaccard | jaccard_max | overlap
    """
    n1, n2 = len(geneids1), len(geneids2)
    n = len(geneids1 & geneids2)
    N = len(nodes)
    if method == "hyper":
	val = stat_utilities.hypergeometric_test_numeric(n, n1, N, n2)
    elif method == "fishers":
	oddsratio, val = stat_utilities.fisher_exact(n, n1 - n, n2 -n, N - n1 - n2 + n, alternative="two-sided")
    elif method == "jaccard":
	val = stat_utilities.jaccard(geneids1, geneids2)
    elif method == "jaccard_max":
	val = stat_utilities.jaccard_max(geneids1, geneids2)
    elif method == "overlap":
	val = n
    else:
	raise ValueError("Uknown method: %s" % method)
    return n, n1, n2, val


##### Proximity related #####

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    #distance = "closest"
    #lengths = network_utilities.get_shortest_path_lengths(network, "../data/toy.sif.pcl")
    #d = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network 
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
	return None # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
	bins = network_utilities.get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
	nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
	nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
	nodes_from, nodes_to = values_random
	#values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
	values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    #pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
	z = 0.0
    else:
	z = (d - m) / s
    return d, z, (m, s) #(z, pval)


def calculate_proximity_multiple(network, from_file=None, to_file=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None, out_file="output.txt"):
    """
    Run proximity on each entries of from and to files in a pairwise manner
    output is saved in out_file (e.g., output.txt)
    """
    nodes = set(network.nodes())
    drug_to_targets, drug_to_category = get_diseasome_genes(from_file, nodes = nodes)
    #drug_to_targets = dict((drug, nodes & targets) for drug, targets in drug_to_targets.iteritems())
    disease_to_genes, disease_to_category = get_diseasome_genes(to_file, nodes = nodes)
    # Calculate proximity values
    print len(drug_to_targets), len(disease_to_genes)
    # Get degree binning 
    bins = network_utilities.get_degree_binning(network, min_bin_size)
    f = open(out_file, 'w')
    f.write("source\ttarget\tn.source\tn.target\td\tz\n")
    for drug, nodes_from in drug_to_targets.iteritems():
	values = []
	for disease, nodes_to in disease_to_genes.iteritems():
	    print drug, disease
	    d, z, (m, s) = calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=bins, n_random=n_random, min_bin_size=min_bin_size, seed=seed, lengths=lengths)
	    values.append((drug, disease, z, len(nodes_from), len(nodes_to), d, m, s))
	    #f.write("%s\t%s\t%f\t%f\t%f\t%f\n" % (drug, disease, z, d, m, s))
	values.sort(key=lambda x: x[2])
	for drug, disease, z, k, l, d, m, s in values:
	    #f.write("%s\t%s\t%f\t%d\t%d\t%f\t%f\t%f\n" % (drug, disease, z, k, l, d, m, s))
	    f.write("%s\t%s\t%d\t%d\t%f\t%f\n" % (drug, disease, k, l, d, z))
    f.close()
    return 


def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
	for node_from in nodes_from:
	    values = []
	    for node_to in nodes_to:
		val = network_utilities.get_shortest_path_length_between(network, node_from, node_to)
		values.append(val)
	    d = min(values)
	    #print d,
	    values_outer.append(d)
    else:
	for node_from in nodes_from:
	    values = []
	    vals = lengths[node_from]
	    for node_to in nodes_to:
		val = vals[node_to]
		values.append(val)
	    d = min(values)
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


### Separation related

def calculate_separation_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    nodes_network = set(network.nodes())
    if len(set(nodes_from) & nodes_network) == 0 or len(set(nodes_to) & nodes_network) == 0:
	return None # At least one of the node group not in network
    d = get_separation(network, nodes_from, nodes_to, lengths)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
	bins = network_utilities.get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
	nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
	nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
	nodes_from, nodes_to = values_random
	values[i] = get_separation(network, nodes_from, nodes_to, lengths)
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
	z = 0.0
    else:
	z = (d - m) / s
    return d, z, (m, s) #(z, pval)


def get_separation(network, nodes_from, nodes_to, lengths=None):
    dAA = numpy.mean(get_separation_within_set(network, nodes_from))
    dBB = numpy.mean(get_separation_within_set(network, nodes_to))
    dAB = numpy.mean(get_separation_between_sets(network, nodes_from, nodes_to))
    d = dAB - (dAA + dBB) / 2.0
    return d


def get_separation_between_sets(network, nodes_from, nodes_to, lengths=None):
    """
    Calculate dAB in separation metric proposed by Menche et al. 2015
    """
    values = []
    target_to_values = {}
    source_to_values = {}
    for source_id in nodes_from:
	for target_id in nodes_to:
	    if lengths is not None:
		d = lengths[source_id][target_id] 
	    else:
		d = network_utilities.get_shortest_path_length_between(network, source_id, target_id)
	    source_to_values.setdefault(source_id, []).append(d)
	    target_to_values.setdefault(target_id, []).append(d)
    # Distances to closest node in nodes_to (B) from nodes_from (A)
    for source_id in nodes_from:
	inner_values = source_to_values[source_id]
	values.append(numpy.min(inner_values))
    # Distances to closest node in nodes_from (A) from nodes_to (B)
    for target_id in nodes_to:
	inner_values = target_to_values[target_id]
	values.append(numpy.min(inner_values))
    return values


def get_separation_within_set(network, nodes_from, lengths=None):
    """
    Calculate dAA or dBB in separation metric proposed by Menche et al. 2015
    """
    if len(nodes_from) == 1:
	return [ 0 ]
    values = []
    # Distance to closest node within the set (A or B)
    for source_id in nodes_from:
	inner_values = []
	for target_id in nodes_from:
	    if source_id == target_id:
		continue
	    if lengths is not None:
		d = lengths[source_id][target_id] 
	    else:
		d = network_utilities.get_shortest_path_length_between(network, source_id, target_id)
	    inner_values.append(d)
	values.append(numpy.min(inner_values))
    return values



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


def run_guild(phenotype, node_to_score, network_nodes, network_file, output_dir, executable_path = None, background_score = 0.01, qname=None, method='s'): 
    # Create node file
    node_file = "%s%s.node" % (output_dir, phenotype)
    create_node_file(node_to_score, network_nodes, node_file, background_score)
    output_file = "%s%s.n%s" % (output_dir, phenotype, method)
    # Get and run the GUILD command
    #print strftime("%H:%M:%S - %d %b %Y") #, score_command
    if method == 's': 
	n_repetition = 3 
	n_iteration = 2 
	score_command = ' -s s -n "%s" -e "%s" -o "%s" -r %d -i %d' % (node_file, network_file, output_file, n_repetition, n_iteration)
    elif method == 'd':
	score_command = ' -s d -n "%s" -e "%s" -o "%s"' % (node_file, network_file, output_file)
    elif method == 'r':
	n_iteration = 50 
	score_command = ' -s r -n "%s" -e "%s" -o "%s" -i %d' % (node_file, network_file, output_file, n_iteration)
    elif method == 'p':
	score_command = ' "%s" "%s" "%s" 1' % (node_file, network_file, output_file)
    elif method == 'w':
	score_command = ' "%s" "%s" "%s"' % (node_file, network_file, output_file)
    else:
	raise NotImplementedError("method %s" % method) 
    if qname is None:
	if executable_path is None:
	    if method in ["s", "r", "d"]:
		executable_path = "guild" # assuming accessible guild executable
	    else:
		executable_path = "netwalk.sh" # assuming R and netwalk.sh is accessible  
	score_command = executable_path + score_command
	print score_command
	os.system(score_command)
    else:
	#os.system("qsub -cwd -o out -e err -q %s -N %s -b y %s" % (qname, scoring_type, score_command))
	#print "qsub -cwd -o out -e err -q %s -N guild_%s -b y %s" % (qname, drug, score_command)
	print "%s" % (score_command.replace('"', ''))
    return score_command


def guildify_multiple(network_file, to_file, output_dir, from_file=None, out_file="guild.txt", method="s", executable_path=None):
    """
    to_file: seeds
    If from_file is not None, returns a dictionary containing average z scores of targets to source, otherwise returns empty dictionary
    method: d | s | r | w | p
    (netshort | netscore | page rank | random walk | propagation)
    """
    if from_file is not None and os.path.exists(out_file):
	target_to_source_score = dict(line.strip("\n").split() for line in open(out_file).readlines())
	return target_to_source_score 
    target_to_source_score = {}
    network = get_network(network_file, only_lcc = True) # using LCC 
    if network_file.endswith(".lcc"):
	network_lcc_file = network_file
    else:
	network_lcc_file = network_file + ".lcc"
    nodes = set(network.nodes())
    disease_to_genes, disease_to_category = get_diseasome_genes(to_file, nodes = nodes)
    if not os.path.exists(output_dir):
	print "Creating output directory", output_dir
	os.makedirs(output_dir)
    # Generate background file (for P-value calculation)
    if not os.path.exists(output_dir + "/background.node"):
	node_to_degree = dict(network.degree())
	n = max(map(len, disease_to_genes.values()))
	values = node_to_degree.items()
	values.sort(key=lambda x: -x[1])
	#k = 1.0 * max(node_to_degree.values())
	values = set(zip(*values[:n])[0])
	f = open(output_dir + "/background.node", 'w')
	for node, degree in node_to_degree.iteritems():
	    #score = degree/k
	    if node in values: score = 1
	    else: score = 0.01
	    f.write("%s %f\n" % (node, score))
	f.close()
    if from_file is not None:
	drug_to_targets, drug_to_category = get_diseasome_genes(from_file, nodes = nodes)
	f = open(out_file, 'w')
	f.write("source\ttarget\tscore\n")
    for target, geneids in disease_to_genes.iteritems():
	#print target, len(geneids)
	target_mod = text_utilities.convert_to_R_string(target)
	target_to_score = dict((gene, 1.0) for gene in geneids)
	node_file = output_dir + "%s.n%s" % (target_mod, method) 
	if os.path.exists(node_file):
	    print "Skipping existing:", node_file
	    continue
	run_guild(target_mod, target_to_score, nodes, network_lcc_file, output_dir, executable_path, background_score = 0.01, qname = "print", method = method) #!
	node_to_score = dict(line.strip("\n").split() for line in open(node_file).readlines())
	if from_file is not None:
	    values = map(float, numpy.array(node_to_score.values()))
	    m = numpy.mean(values)
	    s = numpy.std(values)
	    for source, geneids in drug_to_targets.iteritems():
		score = -numpy.mean([(float(node_to_score[gene]) - m) / s for gene in geneids])
		f.write("%s\t%s\t%f\n" % (source, target, score))
		d = target_to_source_score.setdefault(target, {})
		d[source] = score
    if from_file is not None:
	f.close()
    return target_to_source_score


def get_scores(score_file):
    """
	Parses scores from a scoring file created by GUILD (node <whitespace> score), returns a dictionary where the values are floats.
    """
    nodes, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = score_file, store_edge_type = False, delim=None, data_to_float=True)
    return node_to_score



### DIAMOnD related ###

def get_diamond_genes(network_file, seeds, file_name, only_lcc=True):
    network = get_network(network_file, only_lcc=only_lcc) 
    nodes = set(network.nodes())
    seeds = set(seeds) & nodes
    #print len(seeds)
    n_iteration = 500
    if not os.path.exists(file_name):
	diamond.DIAMOnD(network, seeds, n_iteration, alpha = 1, outfile = file_name)
    f = open(file_name)
    f.readline()
    genes = []
    for line in f:
	rank, geneid = line.strip("\n").split()
	genes.append(geneid)
    f.close()
    if not os.path.exists(file_name + ".coverage"):
	f_out = open(file_name + ".coverage", 'w')
	n = float(len(seeds))
	component = network.subgraph(seeds)
	#component = max(networkx.connected_components(component), key=len)
	components = max(network_utilities.get_connected_components(network, False), key=len)
	f_out.write("%s %f\n" % ("0", len(component & seeds)/n))
	for i, gene in enumerate(genes):
	    rank = i + 1
	    component = network.subgraph(genes[:rank] + list(seeds))
	    #component = max(networkx.connected_components(component), key=len)
	    components = max(network_utilities.get_connected_components(network, False), key=len)
	    f_out.write("%s %f\n" % (rank, len(component & seeds)/n))
	f_out.close()
    return genes, nodes


### Functional enrichment related ###

def check_functional_enrichment(id_list, background_id_weights = None, id_type = "genesymbol", species = "Homo sapiens", mode="unordered", evidences = None, out_file_name = None, tex_format = False):
    """
    id_type = "geneid" # "uniprotacession" # "genesymbol"
    evidences = ['EXP', 'IDA', 'IEP', 'IGI', 'IMP', 'ISA', 'ISM', 'ISO', 'ISS', 'IGC'] # 'IPI'
    evidences = None corresponds to ['EXP', 'IC', 'IDA', 'IEA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NAS', 'RCA', 'TAS']
    for custom associations: association = [["GO:0006509", "351"], ["GO:0048167", "348", "5663", "5664", "23621"], ["GO:0097458", "1005", "1006", "1007"], ["GO:0048487", "1", "2", "351"], ["GO:0048488", map(str, range(1000,2000))]]
    for backgroud id weights, such as occurrence frequency: (gene, weight) 
    """
    if out_file_name is not None:
	f_output = open(out_file_name, 'w').write
    else:
	from sys import stdout 
	f_output = stdout.write
    return functional_enrichment.check_functional_enrichment(id_list, background_id_weights, id_type, f_output, species = species, mode = mode, tex_format = tex_format, support = evidences, associations = None) 



