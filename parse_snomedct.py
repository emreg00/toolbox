import networkx

def main():
    base_dir = "/home/emre/data/ontology/snomedct/SnomedCT_Release_INT_20130731/RF2Release/Full/Terminology/"
    desc_file = base_dir + "sct2_Description_Full-en_INT_20130731.txt"
    rel_file = base_dir + "sct2_Relationship_Full_INT_20130731.txt"
    s = SNOMEDCT(desc_file, rel_file)
    #concept = "Pharmaceutical / biologic product"
    #print s.get_concept_ids(concept)
    g = s.get_ontology()
    #print len(g.nodes()), len(g.edges())
    concept = "Triazole antifungals"
    concept_id = s.get_concept_ids(concept)[1]
    print concept_id, g.edges([concept_id])
    for u, v in g.edges([concept_id]):
	print v, s.get_concept(v)
    return

class SNOMEDCT(object):

    def __init__(self, file_name_desc, file_name_rel):
	self.file_name_desc = file_name_desc
	self.file_name_rel = file_name_rel
	self.delim = "\t"
	self.ontology = None 
	return

    def get_concept_ids(self, concept):
	f = open(self.file_name_desc)
	header = f.readline().strip("\n")
	col_to_idx = dict((val.lower(), i) for i, val in enumerate(header.split(self.delim)))
	concept_ids = None
	for line in f:
	    words = line.strip("\n").split(self.delim)
	    if words[col_to_idx["term"]] == concept:
		concept_id = words[col_to_idx["conceptid"]]
		if concept_ids is None:
		    concept_ids = [concept_id]
		else:
		    concept_ids.append(concept_id)
	if concept_ids is None:
	    raise ValueError("Concept not found")
	return concept_ids

    def get_concept(self, concept_id):
	f = open(self.file_name_desc)
	header = f.readline().strip("\n")
	col_to_idx = dict((val.lower(), i) for i, val in enumerate(header.split(self.delim)))
	concept = None
	for line in f:
	    words = line.strip("\n").split(self.delim)
	    if words[col_to_idx["conceptid"]] == concept_id:
		concept = words[col_to_idx["term"]]
		break
	if concept is None:
	    raise ValueError("Concept not found")
	return concept

    def get_ontology(self, root_concept="Pharmaceutical / biologic product", relation_types=["Is a"]):
	if self.ontology is not None:
	    return self.ontology
	self.ontology = networkx.DiGraph()
	root = self.get_concept_ids(root_concept)[0]
	valid_relations = set([self.get_concept_ids(relation)[0] for relation in relation_types])
	f = open(self.file_name_rel)
	header = f.readline().strip("\n")
	col_to_idx = dict((val.lower(), i) for i, val in enumerate(header.split(self.delim)))
	for line in f:
	    words = line.strip("\n").split(self.delim)
	    if words[col_to_idx["typeid"]] in valid_relations:
		source_id = words[col_to_idx["sourceid"]]
		target_id = words[col_to_idx["destinationid"]]
		self.ontology.add_edge(target_id, source_id)
	self.ontology = networkx.dfs_tree(self.ontology, root)
	return self.ontology


if __name__ == "__main__":
    main()

