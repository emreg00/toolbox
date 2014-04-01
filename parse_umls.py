import networkx

def main():
    base_dir = "/home/emre/arastirma/data/umls/2013AA/META/"
    desc_file = base_dir + "MRCONSO.RRF"
    rel_file = base_dir + "MRREL.RRF"
    #get_basic_info(desc_file, rel_file)
    get_drug_info(desc_file, rel_file)
    return

class UMLS(object):

    def __init__(self, file_name_desc, file_name_rel):
	self.file_name_desc = file_name_desc
	self.file_name_rel = file_name_rel
	self.delim = "|"
	self.ontology = None 
	self.concept_id_to_values = None
	self.concept_to_concept_ids = None
	self.concept_id_to_relations = None
	return

    def get_concept_info(self, concept_types=None, concept_sources=None):
	if self.concept_id_to_values is not None and self.concept_to_concept_id is not None:
	    return self.concept_id_to_values, self.concept_to_concept_id
	self.concept_id_to_values = {}
	self.concept_to_concept_id = {}
	f = open(self.file_name_desc)
	header_names = ["CUI", "LAT", "TS", "LUI", "STT", "SUI", "ISPREF", "AUI", "SAUI", "SCUI", "SDUI", "SAB", "TTY", "CODE", "STR", "SRL", "SUPPRESS", "CVF", "dummy"]
	# CUI / LAT (ENG) / TS (P) / STT (PF/VO) all / ISPREF (Y) / SCUI - source based id / SAB - source / TTY (PT/SY) pt-prefered sy-synonym / CODE similar to SCUI / STR
	col_to_idx = dict((val.lower(), i) for i, val in enumerate(header_names))
	for line in f:
	    words = line.strip("\n").split(self.delim)
	    concept_id = words[col_to_idx["cui"]]
	    #if concept_id == "C0360380":
	    #	print len(words), words
	    #	print words[col_to_idx["ts"]], words[col_to_idx["ispref"]], words[col_to_idx["tty"]]
	    if words[col_to_idx["lat"]] != "ENG" or words[col_to_idx["ispref"]] != "Y": # words[col_to_idx["ts"]] != "P" 
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
	concept_to_values, concept_to_concept_id = self.get_concept_info()
	return concept_to_concept_id[concept]

    def get_concepts(self, concept_id, concept_types=set(["PT", "PN", "EN", "EP", "FN", "SY", "PM"])):
	concepts = []
	concept_to_values, concept_to_concept_id = self.get_concept_info()
	values = concept_to_values[concept_id]
	for source, vals in values.iteritems():
	    for concept, source_id, concept_type in vals:
		if concept_type in concept_types:
		    concepts.append((source, concept, concept_type))
		#else:
		#    print concept_type
	#if len(concepts) == 0:
	#    raise ValueError("Concept not found")
	return concepts

    def get_relations(self, relation_types=set(["treats", "treated_by"]), source_types=None): # , "may_treat", "may_be_treated"
	if self.concept_id_to_relations is not None:
	    return self.concept_id_to_relations 
	self.concept_id_to_relations = {}
	f = open(self.file_name_rel)
	header_names = ["CUI1", "AUI1", "STYPE1", "REL", "CUI2", "AUI2", "STYPE2", "RELA", "RUI", "SRUI", "SAB", "SL", "RG", "DIR", "SUPPRESS", "CVF", "dummy"]
	col_to_idx = dict((val.lower(), i) for i, val in enumerate(header_names))
	for line in f:
	    words = line.strip("\n").split(self.delim)
	    relation = words[col_to_idx["rela"]]
	    if relation not in relation_types:
		continue
	    source_id = words[col_to_idx["cui1"]]
	    target_id = words[col_to_idx["cui2"]]
	    source = words[col_to_idx["sab"]]
	    if source_types is not None and source not in source_types: 
		continue
	    d = self.concept_id_to_relations.setdefault(target_id, {})
	    d[source_id] = (relation, source)
	return self.concept_id_to_relations 

    def get_ontology(self, root_concept="Pharmaceutical / biologic product", relation_types=set(["isa"]), source_types=set(["SNOMEDCT"])):
	if self.ontology is not None:
	    return self.ontology
	self.ontology = networkx.DiGraph()
	root = self.get_concept_id(root_concept)
	f = open(self.file_name_rel)
	header_names = ["CUI1", "AUI1", "STYPE1", "REL", "CUI2", "AUI2", "STYPE2", "RELA", "RUI", "SRUI", "SAB", "SL", "RG", "DIR", "SUPPRESS", "CVF", "dummy"]
	col_to_idx = dict((val.lower(), i) for i, val in enumerate(header_names))
	i = 0
	for line in f:
	    words = line.strip("\n").split(self.delim)
	    source_id = words[col_to_idx["cui1"]]
	    target_id = words[col_to_idx["cui2"]]
	    relation = words[col_to_idx["rela"]]
	    source = words[col_to_idx["sab"]]
	    #if source_id == root or target_id == root:
	    #	print self.get_concepts(source_id), relation, self.get_concepts(target_id), source
	    if relation not in relation_types:
		continue
	    if source_types is not None and source not in source_types: 
		continue
	    self.ontology.add_edge(target_id, source_id)
	    i += 1
	    #if i > 1000:
	    #	break
	#print len(self.ontology.nodes()), len(self.ontology.edges())
	networkx.dfs_tree(self.ontology.reverse(), root)
	return self.ontology

    def get_drug_disease_relations(self):
	drug_to_diseases = {}
	for nodes in self.get_ontology().edges():
	    for node in nodes:
		try:
		    rels = self.get_relations()[node]
		except:
		    continue
		for cid, values in rels.iteritems():
		    relation, source = values
		    #if relation != "treats":
		    #    continue
		    for source, concept, concept_type in self.get_concepts(node):
			for source2, concept2, concept_type2 in self.get_concepts(cid):
			    drug_to_diseases.setdefault(concept, set()).add(concept2)
	return drug_to_diseases


def get_mesh_id_to_name(desc_file):
    u = UMLS(desc_file, None)
    concept_id_to_values, concept_to_concept_id = u.get_concept_info(concept_types = set(["MH"]), concept_sources = set(["MSH"]))
    source_id_to_concept = {}
    for concept_id, values in concept_id_to_values.iteritems():
    	for concept, source_id, concept_type in values["MSH"]:
	    source_id_to_concept[source_id] = concept
    return source_id_to_concept

def get_basic_info(desc_file, rel_file):
    u = UMLS(desc_file, rel_file)
    concept = "Triazole antifungals"
    concept_id = u.get_concept_id(concept)
    print concept_id
    concepts = u.get_concepts(concept_id)
    print concepts
    g = u.get_ontology()
    print len(g.nodes()), len(g.edges())
    print concept_id, g.edges([concept_id])
    for s, v in g.edges([concept_id]):
	print u.get_concepts(v)
    concept_id = "C0360363" # Azole antifungal
    rels = u.get_relations()[concept_id]
    for cid, values in rels.iteritems():
	print cid, values
    return

def get_drug_info(desc_file, rel_file):
    u = UMLS(desc_file, rel_file)
    drug_to_diseases = u.get_drug_disease_relations()
    for drug, diseases in drug_to_diseases.iteritems():
    	print drug, diseases
    return

def get_disease_specific_drugs(umls, name_to_drug, synonym_to_drug, phenotypes):
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

