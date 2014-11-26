import networkx, os, cPickle

def main():
    base_dir = "/home/emre/arastirma/data/ontology/umls/2013AA/META/"
    desc_file = base_dir + "MRCONSO.RRF"
    rel_file = base_dir + "MRREL.RRF"
    #g = get_mesh_disease_ontology(desc_file, rel_file)
    #get_basic_info(desc_file, rel_file)
    #get_drug_info(desc_file, rel_file)
    mesh_id_to_name, concept_id_to_mesh_id, mesh_id_to_name_with_synonyms = get_mesh_id_mapping(desc_file, rel_file)
    print "mesh dict:", len(mesh_id_to_name), len(mesh_id_to_name_with_synonyms)
    print mesh_id_to_name["D003924"]
    print concept_id_to_mesh_id["C0011860"]
    print mesh_id_to_name_with_synonyms["D003924"]
    for mesh_id in [ "D003924", "D001769", "D005947", "D004493", "D006943" ]:
	print mesh_id, mesh_id in mesh_id_to_name
    return


class UMLS(object):


    def __init__(self, file_name_desc, file_name_rel, concept_types = None, concept_sources = None, only_preferred = False):
	self.file_name_desc = file_name_desc
	self.file_name_rel = file_name_rel
	self.delim = "|"
	self.ontology = None 
	self.concept_id_to_values = None
	self.concept_to_concept_id = None
	self.concept_id_to_relations = None
	self._get_concept_info(concept_types, concept_sources, only_preferred)
	return


    def _get_concept_info(self, concept_types = None, concept_sources = None, only_preferred = False):
	"""
	Parses MRCONSO file to get concept info, typically called without any arguments and saved to the dictionary
	"""
	if self.concept_id_to_values is None and self.concept_to_concept_id is None:
	    self.concept_id_to_values = {}
	    self.concept_to_concept_id = {}
	    f = open(self.file_name_desc)
	    header_names = ["CUI", "LAT", "TS", "LUI", "STT", "SUI", "ISPREF", "AUI", "SAUI", "SCUI", "SDUI", "SAB", "TTY", "CODE", "STR", "SRL", "SUPPRESS", "CVF", "dummy"]
	    # CUI / LAT (ENG) / TS (P) / STT (PF/VO) all / ISPREF (Y) / SCUI - source based id / SAB - source / TTY (PT/SY) pt-preferred sy-synonym / CODE similar to SCUI / STR
	    col_to_idx = dict((val.lower(), i) for i, val in enumerate(header_names))
	    for line in f:
		words = line.strip("\n").split(self.delim)
		concept_id = words[col_to_idx["cui"]]
		#if concept_id == "C0360380":
		#	print len(words), words
		#	print words[col_to_idx["ts"]], words[col_to_idx["ispref"]], words[col_to_idx["tty"]]
		if words[col_to_idx["lat"]] != "ENG": # words[col_to_idx["ts"]] != "P" 
		    continue
		if only_preferred and words[col_to_idx["ispref"]] != "Y": 
		    continue
		concept_type = words[col_to_idx["tty"]]
		if concept_types is not None and concept_type not in concept_types: 
		    continue
		source = words[col_to_idx["sab"]]
		if concept_sources is not None and source not in concept_sources: 
		    continue
		concept = words[col_to_idx["str"]]
		source_id = words[col_to_idx["code"]]
		d = self.concept_id_to_values.setdefault(concept_id, {})
		d.setdefault(source, set()).add((concept, source_id, concept_type))
		if concept_id in self.concept_to_concept_id:
		    print "Concept id conflict - overwriting:", concept, self.concept_to_concept_id[concept], concept_id
		self.concept_to_concept_id[concept] = concept_id
	return self.concept_id_to_values, self.concept_to_concept_id


    def get_concept_id(self, concept):
	return self.concept_to_concept_id[concept]


    def get_values_by_concept_id(self, concept_id):
	return self.concept_id_to_values[concept_id]


    def get_concepts(self, concept_id, concept_sources = None, concept_types = None): 
	concepts = []
	values = self.get_values_by_concept_id(concept_id)
	for source, vals in values.iteritems():
	    if concept_sources is not None and source not in concept_sources:
		continue
	    for concept, source_id, concept_type in vals:
		if concept_types is not None and concept_type not in concept_types:
		    continue
		concepts.append((source, concept, concept_type))
		#else:
		#    print concept_type
	#if len(concepts) == 0:
	#    raise ValueError("Concept not found")
	return concepts


    def get_relations(self, relation_types = None, relation_a_types = None, source_types = None): # , "may_treat", "may_be_treated"
	"""
	Parses MRREL file to get relation info, typically called with relation type parameters and saved to the dictionary
	"""
	if self.concept_id_to_relations is None:
	    self.concept_id_to_relations = {}
	    f = open(self.file_name_rel)
	    header_names = ["CUI1", "AUI1", "STYPE1", "REL", "CUI2", "AUI2", "STYPE2", "RELA", "RUI", "SRUI", "SAB", "SL", "RG", "DIR", "SUPPRESS", "CVF", "dummy"]
	    col_to_idx = dict((val.lower(), i) for i, val in enumerate(header_names))
	    for line in f:
		words = line.strip("\n").split(self.delim)
		relation = words[col_to_idx["rel"]]
		relation_a = words[col_to_idx["rela"]]
		if relation_types is not None and relation not in relation_types:
		    continue
		if relation_a_types is not None and relation_a not in relation_a_types:
		    continue
		source_id = words[col_to_idx["cui1"]]
		target_id = words[col_to_idx["cui2"]]
		source = words[col_to_idx["sab"]]
		if source_types is not None and source not in source_types: 
		    continue
		d = self.concept_id_to_relations.setdefault(target_id, {})
		d[source_id] = (relation, source)
		if relation_a != "":
		    d[source_id] = (relation_a, source)
	return self.concept_id_to_relations 


    def get_ontology(self, root_concept = None, relation_types = None, relation_a_types = None, source_types = None):
	"""
	Gets the graph (ontology tree) from MRREL file, typically called with relation type parameters and saved as a Networkx Graph object
	"""
	if self.ontology is None:
	    self.ontology = networkx.DiGraph()
	    f = open(self.file_name_rel)
	    header_names = ["CUI1", "AUI1", "STYPE1", "REL", "CUI2", "AUI2", "STYPE2", "RELA", "RUI", "SRUI", "SAB", "SL", "RG", "DIR", "SUPPRESS", "CVF", "dummy"]
	    col_to_idx = dict((val.lower(), i) for i, val in enumerate(header_names))
	    i = 0
	    for line in f:
		words = line.strip("\n").split(self.delim)
		source_id = words[col_to_idx["cui1"]]
		target_id = words[col_to_idx["cui2"]]
		relation = words[col_to_idx["rel"]]
		relation_a = words[col_to_idx["rela"]]
		source = words[col_to_idx["sab"]]
		if relation_types is not None and relation not in relation_types:
		    continue
		if relation_a_types is not None and relation_a not in relation_a_types:
		    continue
		if source_types is not None and source not in source_types: 
		    continue
		#if source_id == root or target_id == root:
		#	print self.get_concepts(source_id), relation, self.get_concepts(target_id), source
		self.ontology.add_edge(target_id, source_id)
		i += 1
		#if i > 1000:
		#	break
	    self.ontology.reverse()
	if root_concept is not None:
	    root = self.get_concept_id(root_concept)
	    g = networkx.dfs_tree(self.ontology, root)
	else:
	    g = self.ontology
	return g


    def get_drug_disease_relations(self):
	drug_to_diseases = {}
	concept_types = set(["MH", "PF", "PT", "PN", "EN", "EP", "FN", "SY", "PM"])
	for nodes in self.get_ontology(root_concept = "Pharmaceutical / biologic product", relation_a_types = set(["isa"]), source_types = set(["SNOMEDCT"])).edges():
	    for node in nodes:
		try:
		    rels = self.get_relations(relation_a_types=set(["treats", "may_treat"]), source_types = None)[node]
		except:
		    continue
		for cid, values in rels.iteritems():
		    relation, source = values
		    #if relation != "treats":
		    #    continue
		    for source, concept, concept_type in self.get_concepts(node, concept_types = concept_types):
			for source2, concept2, concept_type2 in self.get_concepts(cid, concept_types = concept_types):
			    drug_to_diseases.setdefault(concept, set()).add(concept2)
	return drug_to_diseases


def get_mesh_id_mapping(desc_file, rel_file, only_diseases = True, dump_file = None): 
    if dump_file is not None and os.path.exists(dump_file):
	values = cPickle.load(open(dump_file))
	source_id_to_concept, concept_id_to_mesh_id, source_id_to_concepts = values
	return source_id_to_concept, concept_id_to_mesh_id, source_id_to_concepts
    umls = UMLS(desc_file, rel_file)
    concept_ids_disease = None
    if only_diseases:
	g = get_mesh_disease_ontology(desc_file, rel_file, umls)
	concept_ids_disease = set(g.nodes())
    source_id_to_concept = {} # only main headers
    source_id_to_concepts = {} # all concepts including synonyms
    concept_id_to_mesh_id = {}
    for concept_id, values in umls.concept_id_to_values.iteritems():
	if concept_ids_disease is not None and concept_id not in concept_ids_disease:
	    continue
    	for concept, source_id, concept_type in values["MSH"]:
	    if concept_type == "MH": # main heading
                source_id_to_concept[source_id] = concept
	    source_id_to_concepts.setdefault(source_id, set()).add(concept)
	    #if concept_id in concept_id_to_mesh_id and concept_id_to_mesh_id[concept_id] != source_id:
	    #	print "Inconsistency", concept_id, source_id 
	    concept_id_to_mesh_id[concept_id] = source_id
    if dump_file is not None:
	values = (source_id_to_concept, concept_id_to_mesh_id, source_id_to_concepts)
	cPickle.dump(values, open(dump_file, 'w'))
    return source_id_to_concept, concept_id_to_mesh_id, source_id_to_concepts


def get_mesh_disease_ontology(desc_file, rel_file, umls = None):
    if umls is None:
	umls = UMLS(desc_file, rel_file)
    root = "Diseases (MeSH Category)"
    sources = set(["MSH"]) 
    relations = set(["PAR"]) 
    g = umls.get_ontology(root_concept = root, relation_types = relations, source_types = sources)
    #print "Shrunk ontology:", len(g.nodes()), len(g.edges())
    #for node in g.neighbors(umls.get_concept_id(root)):
    #	print node, umls.get_concepts(node, concept_sources = sources)
    return g


def get_snomedct_drug_ontology(desc_file, rel_file, umls = None):
    if umls is None:
	umls = UMLS(desc_file, rel_file)
    root = "Pharmaceutical / biologic product"
    sources = set(["SNOMEDCT"])
    relations = set(["isa"])
    g = umls.get_ontology(root_concept = root, relation_types = relations, source_types = sources)
    return g


def get_basic_info(desc_file, rel_file):
    u = UMLS(desc_file, rel_file)
    concept = "Diabetes Mellitus" #"Triazole antifungals"
    sources = set(["MSH"]) # set(["SNOMEDCT"])
    relations = set(["PAR"]) # set(["isa"])
    concept_id = u.get_concept_id(concept)
    print concept, concept_id
    concepts = u.get_concepts(concept_id, concept_sources = sources)
    print concepts
    root = "Diseases (MeSH Category)" #"Pharmaceutical / biologic product"
    g = u.get_ontology(root_concept = root, relation_types = relations, source_types = sources)
    print len(g.nodes()), len(g.edges())
    print concept_id, g.edges([concept_id])
    for s, v in g.edges([concept_id]):
	print s, v, u.get_concepts(v, concept_sources = sources)
    concept_id = "C0011849" #"C0360363" # Azole antifungal
    rels = u.get_relations(relation_types = relations, source_types = sources)[concept_id] 
    for cid, values in rels.iteritems():
	print cid, values
    return


def get_drug_info(desc_file, rel_file):
    u = UMLS(desc_file, rel_file)
    drug_to_diseases = u.get_drug_disease_relations()
    for drug, diseases in drug_to_diseases.iteritems():
    	print drug, diseases
    return


def get_disease_specific_drugs(umls, selected_drugs, name_to_drug, synonym_to_drug, phenotypes):
    drug_to_diseases = umls.get_drug_disease_relations()
    disease_to_drugs = {}
    for drug, diseases in drug_to_diseases.iteritems():
	drug = drug.split()[0].lower() 
	if drug in name_to_drug:
	    drugbank_id = name_to_drug[drug]
	elif drug in synonym_to_drug:
	    drugbank_id = synonym_to_drug[drug]
	else:
	    continue
	if drugbank_id not in selected_drugs:
	    continue
	for description in diseases:
	    description = description.lower()
	    for phenotype in phenotypes:
		disease_mod = phenotype.replace(" and ", ", ")
		phrases = disease_mod.split(",")
		values = []
		for phrase in phrases:
		    inner_values = []
		    words = phrase.strip().split()
		    for i, token in enumerate(words):
			if token.endswith("'s"):
			    token = token[:-2]
			if i == len(words) - 1:
			    if token[-1] == "s":
				token = token[:-1]
			if token in ("disease", "disorder", "syndrome"):
			    continue
			inner_values.append(token)
		    #if len(inner_values) > 0:
		    values.append(" ".join(inner_values))
		if all([ description.find(word.strip()) != -1 for word in values ]): # phenotype.split(",")
		    disease_to_drugs.setdefault(phenotype, set()).add(drugbank_id)
    return disease_to_drugs


def old_get_disease_specific_drugs(umls, drug_to_name, phenotypes):
    import re
    #drug_to_diseases = {"telmisartan": set(['Diabetic renal disease', 'congestive cardiac failure', 'congestive heart failure chf', 'left ventricular dysfunction', 'HBP', 'failure congestive heart']) 
    drug_to_diseases = umls.get_drug_disease_relations()
    exps = [ re.compile(keyword.lower()) for keyword in phenotypes ]
    drug_id_to_exp = {} 
    for drug_id, keyword in drug_to_name.iteritems():
	try:
	    for l in "[{}]":
		keyword = keyword.replace(l, "_")
	    exp = re.compile(keyword.lower())
	except:
	    print keyword
	    continue
	drug_id_to_exp[drug_id] = exp
    disease_to_drugs = {}
    for drug, diseases in drug_to_diseases.iteritems():
	#drug = drug.lower().split()[0]
	#print drug, diseases
	drugbank_ids = [] #None 
	for drug_id, drug_name in drug_to_name.iteritems():
	    if drug_id not in drug_id_to_exp:
		continue
	    exp_drug = drug_id_to_exp[drug_id]
	    if exp_drug.search(drug.lower()) is not None:
		#if len(drugbank_ids) > 0: # is not None:
		    #raise ValueError("Duplicate match for drug " + drug_id)
		    #print "Duplicate match for drug ", drug, drug_id, drugbank_ids
		drugbank_ids.append(drug_id)
	if len(drugbank_ids) == 0:
	    continue
	for disease, exp in zip(phenotypes, exps):
	    if any(map(lambda x: x is not None, [ exp.search(description.lower()) for description in diseases ])):
		selected_drugbank_id = None
		length = 0
		for drugbank_id in drugbank_ids:
		    # choose drug with longer name 
		    val = len(drug_to_name[drugbank_id])
		    if val > length:
			selected_drugbank_id = drugbank_id
			length = val
		#if len(drugbank_ids) > 1:
		#    print selected_drugbank_id, drugbank_ids
		disease_to_drugs.setdefault(disease, set()).add(selected_drugbank_id)
    return disease_to_drugs



if __name__ == "__main__":
    main()

