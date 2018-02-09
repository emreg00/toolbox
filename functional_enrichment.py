try:
    #from external.funcassociate import client
    #from funcassociate import client
    #from toolbox import func_associate as client
    import func_associate as client
except:
    client = None
    print "Import error: Funcassociate. Make sure that funcassociate is in toolbox!"

def main():
    file_name = "gene_list.txt"
    out_file_name = "gene_list_functions.txt"
    f=open(file_name)
    lines = f.readlines()
    a=map(lambda x: x.strip(), lines)
    check_functional_enrichment(a, None, "genesymbol", open(out_file_name, 'w').write, tex_format = False) 
    return

def check_functional_enrichment_of_human_gene_symbols(file_name, out_file_name):
    f=open(file_name)
    a=map(lambda x: x.strip(), f.readlines())
    check_functional_enrichment(a, None, "genesymbol", open(out_file_name, 'w').write, species = "Homo sapiens", mode = "unordered", request_info=False, tex_format = False) 
    return

def check_functional_enrichment(subset_gene_ids, gene_weights, id_type, output_method=None, species = "Homo sapiens", mode = "unordered", request_info=False, tex_format=False, support=None, associations=None):
    """
	Check GO functional enrichment using funcassociate web service
	subset_gene_ids is a list of gene symbols (without whitespace) or gene ids
	id_type: geneid | genesymbol | uniprotaccession | ...
	species: Homo sapiens | Mus musculus | Rattus norvegicus | Saccharomyces cerevisiae | Caenorhabditis elegans | ...
	support types: ['EXP', 'IC', 'IDA', 'IEA', 'IEP', 'IGC', 'IGI', 'IMP', 'IPI', 'ISA', 'ISM', 'ISO', 'ISS', 'NAS', 'RCA', 'TAS']
    """
    reps = 2000
    client_funcassociate = client.FuncassociateClient()
    if id_type == "geneid":
	id_type = "entrezgene"
    elif id_type == "genesymbol":
	if species == "Homo sapiens":
	    id_type = "hgnc_symbol"
	elif species == "Mus musculus":
	    id_type = "mgi_symbol"
	else:
	    print client_funcassociate.available_namespaces(species=[species]) 
	    raise ValueError("Currently human and mouse symbols are supported!")
    elif id_type == "uniprot": # "uniprotaccession"
	#id_type = "uniprot_accession"
	id_type = "uniprot_swissprot"
    #elif id_type == "uniprotentry": 
    #	id_type = "uniprot_id"
    elif id_type == "sgd":
	id_type = "sgd_systematic"
    else:
	raise ValueError("Unrecognized id_type: %s" % id_type)
    response = client_funcassociate.functionate(query = subset_gene_ids,
                             species = species,
                             namespace = id_type,
                             genespace = gene_weights,
                             mode = mode,
                             reps = reps,
			     support = support,
			     associations = associations)

    if output_method is None:
	return response["over"]

    #headers = ["N", "M", "X", "LOD", "P", "P_adj", "attrib ID", "attrib name"]
    headers = [ "# of genes", "# of genes in the query", "# of total genes", "Log of odds ratio", "P-value", "Adjusted p-value", "GO term ID", "Go term name" ]

    #if mode == "unordered":
    #	headers.pop(1)
    headers.pop(1) # Now that column is always present independent of the mode
    if tex_format:
	output_method("%s\\\\\n" % " & ".join(headers))
    else:
	output_method("%s\n" % "\t".join(headers))

    zero = "< %f" % (1.0/float(reps))

    for row in response["over"]:
	if mode == "unordered":
	    row = row[:1] + row[2:] #row.pop(1)
        if row[4] is 0:
            row[4] = zero
	if mode == "unordered":
	    interval = range(2,5)
	else:
	    interval = range(3,6)
	#print row
	for i in interval:
	    if isinstance(row[i], str) and row[i].startswith("<"):
		#print row[i]
		val = float(row[i].lstrip("<"))
		if tex_format:
		    row[i] = "$<$%.5f" % val
		else:
		    row[i] = "<%.5f" % val
	    else:
		row[i] = "%.5f" % row[i]
	if tex_format:
	    output_method("%s\\\\\n" % " & ".join(map(str, row)))
	else:
	    output_method("%s\n" % "\t".join(map(str, row)))

    if request_info:
	output_method("\nREQUEST INFO\n")
	info = response["request_info"]
	for k in info.keys():
	    output_method("%s: %s\n" % (k, info[k]))
    return response["over"]


# Not recommended to be used when two sets of GO terms are going to be compared
# In such cases redundancy can be removed using GoSemSim R package
def remove_parent_terms(go_terms, g):
    to_remove = set()
    while True:
	for go_term in go_terms:
	    parent_terms = g.edges(go_term)
	    #print set(zip(*parent_terms)[1])
	    if len(parent_terms) != 0:
		to_remove |= go_terms & set(zip(*parent_terms)[1])
	    #print to_remove
	if len(to_remove) == 0:
	    break
	#print len(go_terms & to_remove)
	go_terms -= to_remove
	to_remove = set()
    return go_terms

# For GoSemSim calculation in R
def output_go_terms_and_levels(go_terms, go, output_file, root_id="GO:0008150"):
    """
	root_id = "GO:0008150" # BP
    """
    from networkx import bidirectional_shortest_path
    f_out = open(output_file, 'w')
    f_out.write("go level\n")
    for go_id in go_terms:
	level = len(bidirectional_shortest_path(go, go_id, root_id))
	f_out.write("%s %d\n" % (go_id, level)) 
    f_out.close()
    return


def get_go_ontology(file_name):
    from toolbox import OboParser
    go = OboParser.getOboGraph(file_name)
    return go


def get_functional_enrichment(enrichment_file, go, remove_parents=False, only_biological_processes=False, only_slim=False, logodds_cutoff=0):
    """
	Read functional enrichment file.
	If there are multiple functional enrichment analyses it takes the comment as the key and returns
	a dictionary containing name - go_term pairs. If there is only one analysis, returns the go_terms.
    """

    #from toolbox import OboParser
    #g=OboParser.getOboGraph("/home/emre/arastirma/celldiff/data/GO/gene_ontology.1_2.obo")
    g=go

    go_terms = None
    name = None
    name_to_go_terms ={}

    altid_to_goid = {}
    for goid, data in g.nodes(data=True):
	for altid in data["x"]:
	    altid_to_goid[altid] = goid

    f = open(enrichment_file)
    for line in f:
	line = line.strip()
	if line.startswith("# of"):
	    go_terms = set()
	elif line.startswith("#"):
	    if go_terms is not None: 
		if name is None:
		    name = "generic"
		if remove_parents:
		    go_terms = remove_parent_terms(go_terms, g)
		name_to_go_terms[name] = go_terms
		#print name, go_terms
	    name = line 
	else:
	    words = line.split("\t")
	    try:
		n = int(words[0])
		#pval = float(words[4]) # contains text like <0.0067
	    except:
		print words
		continue
	    go_term = words[5]
	    lodds = float(words[2])
	    if lodds < logodds_cutoff:
		continue
	    if go_term in altid_to_goid:
		go_term = altid_to_goid[go_term]
	    if only_biological_processes:
		if g.node[go_term]['t'] == "biological_process": # and 'a' in g.node[go_term]: #is bp and slim # ("molecular_function", "biological_process"):
		    if only_slim:
			if 'a' in g.node[go_term] and g.node[go_term]['a']:
			    go_terms.add(go_term)
		    else:
			go_terms.add(go_term)
	    else:
		if only_slim:
		    if 'a' in g.node[go_term] and g.node[go_term]['a']:
			go_terms.add(go_term)
		else:
		    go_terms.add(go_term)

    f.close()

    if name is None:
	name = "generic"
    if remove_parents:
	go_terms = remove_parent_terms(go_terms, g)
    name_to_go_terms[name] = go_terms
    #print name, go_terms

    if "generic" in name_to_go_terms:
	return name_to_go_terms["generic"]

    return name_to_go_terms


if __name__ == "__main__":
    main()

