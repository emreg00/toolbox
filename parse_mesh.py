import networkx

def main():
    base_dir = "/home/emre/data/ontology/mesh/"
    rel_file = base_dir + "mtrees2013.bin"
    m = MESH(rel_file)
    g = m.get_ontology()
    #print len(g.nodes()), len(g.edges())
    #concept = "Mycoses"
    concept = "Diabetes Mellitus, Type 2"
    concept_ids = m.get_concept_ids(concept)
    for concept_id in concept_ids:
	print concept_id, g.edges([concept_id])
	for u, v in g.edges([concept_id]):
	    print v, m.get_concept(v)
    return

def get_mesh_term_to_tree_ids(rel_file):
    m = MESH(rel_file)
    g = m.get_ontology()
    return m.concept_to_concept_ids

class MESH(object):

    def __init__(self, file_name):
	self.file_name = file_name
	self.delim = ";"
	self.ontology = None 
	self.concept_to_concept_ids = None
	self.concept_id_to_concept = None
	self.root_concept_ids = None
	return

    def get_concept_ids(self, concept):
	return self.concept_to_concept_ids[concept]

    def get_concept(self, concept_id):
	return self.concept_id_to_concept[concept_id]

    def get_ontology(self):
	if self.ontology is not None:
	    return self.ontology
	self.ontology = networkx.DiGraph()
	self.concept_to_concept_ids = {}
	self.concept_id_to_concept = {}
	self.root_concept_ids = []
	f = open(self.file_name)
	for line in f:
	    words = line.strip("\n").split(self.delim)
	    concept = words[0]
	    concept_id = words[1]
	    self.concept_id_to_concept[concept_id] = concept
	    self.concept_to_concept_ids.setdefault(concept, []).append(concept_id)
	    inner_words = concept_id.split(".")
	    if len(inner_words) == 1:
		self.root_concept_ids.append(inner_words[0])
	    else:
		parent_id = ".".join(inner_words[:-1])
		self.ontology.add_edge(parent_id, concept_id)
	return self.ontology


if __name__ == "__main__":
    main()

