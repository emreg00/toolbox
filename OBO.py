
class OBO(object):

    def __init__(self, file_name, save_synonyms = False):
	from Parsers import OboParser
	self.g = OboParser.getOboGraph(file_name, save_synonyms)
	self.precalculated_descendants = {}
	self.child_to_parent = None
	return

    def get_term(self, v):
	return self.g.node[v]

    def get_term_name(self, v):
	return self.g.node[v]['n']

    def get_term_synonyms(self, v):
	return self.g.node[v]['s']

    def get_term_relations(self, v):
	#return self.g.neighbors(v)
	return [ (u, data['r']) for u, data in self.g[v].iteritems() ]

    def get_term_relation_dict(self, v):
	return self.g[v]

    def get_descendants(self, id):
	"""
	Gets all the descendants
	"""
	if self.precalculated_descendants.has_key(id):
	    return self.precalculated_descendants[id]
	result = set()
	for current_descendant_id in self.g.neighbors(id):
	    if current_descendant_id == id:
		return result
	    else:
		if current_descendant_id not in result:
		    result.add(current_descendant_id)
		    result.update(self.get_descendants(current_descendant_id))
	self.precalculated_descendants[id] = result
	return result

    def get_ontology_extended_id_mapping(self, terms=None):
	if terms is not None:
	    nodes = terms
	else:
	    if self.child_to_parent is not None:
		return self.child_to_parent
	    nodes = self.g.nodes()

	for id in nodes:
	    self.get_descendants(id)
	
	self.child_to_parent = {}
	for key, values in self.precalculated_descendants.iteritems():
	    for value in values:
		self.child_to_parent.setdefault(value, set()).add(key)

	return self.child_to_parent

    def get_nested_ontology_mapping(self, from_ontology_prefix,
	    to_ontology_prefix):
	hsdl_to_hge = {} # hsdl_to_go
	for v, data in self.g.nodes(data=True):
	    if not v.startswith(to_ontology_prefix): #"HGE:"):
		continue
	    relations = self.g[v]
	    while len(relations) > 0:
		for child, values in relations.iteritems(): # Currently tracing only the last child
		    if from_ontology_prefix == "*" or child.startswith(from_ontology_prefix): #"HSDL:"):
			#type = values['r']
			hsdl_to_hge.setdefault(child, set()).add(v)
			relations = {}
		    else:
			relations = self.g[child]
	#for key, value in hsdl_to_hge.iteritems():
	#	print key, sorted(list(value))
	return hsdl_to_hge



