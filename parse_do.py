import OBO

def main():
    disease_ontology_file = "../../data/ontology/doid.obo"
    icd_to_mesh_ids = get_icd_to_mesh_ids(disease_ontology_file, id_type="ICD9CM")
    return


def get_icd_to_mesh_ids(disease_ontology_file, id_type="ICD9CM"):
    """
    id_type: ICD9CM | ICD10CM
    """
    name_to_id, id_to_mesh_ids, mesh_id_to_type_to_ids = get_do_mesh_id_mapping(disease_ontology_file)
    #print mesh_id_to_type_to_ids.items()[:5] 
    icd_to_mesh_ids = {}
    for mesh_id, type_to_ids in mesh_id_to_type_to_ids.iteritems():
	if id_type in type_to_ids:
	    for val in type_to_ids[id_type]:
		words = val.split("-")
		if len(words) > 1:
		    icd1 = int(words[0].split(".")[0])
		    icd2 = int(words[1].split(".")[0])
		    icds = map(str, range(icd1, icd2+1))
		else:
		    icds = [ words[0].split(".")[0] ]
		for icd in icds:
		    icd_to_mesh_ids.setdefault(icd, set()).add(mesh_id)
    return icd_to_mesh_ids


def get_do_mesh_id_mapping(do_file):
    do = OBO.OBO(do_file, save_synonyms = True)
    name_to_id = {}
    id_to_mesh_ids = {}
    mesh_id_to_type_to_ids = {}
    for node, data in do.g.nodes_iter(data=True):
	name = data['n']
	name_to_id[name] = node
	if 's' in data:
	    for name in data['s']:
		name_to_id[name] = node
	id_types = ["OMIM", "ICD9CM", "ICD10CM"] #"MESH"
	id_dict = {} #dict((id_type, []) for id_type in id_types) 
	mesh_ids = []
	for xref in data["xref"]:
	    vals = xref.split(":")
	    if len(vals) != 2:
		print "Ontology format inconsistency!", vals
		continue
	    id_type, id_val = vals
	    if id_type == "MESH":
		mesh_ids.append(id_val)
	    elif id_type in id_types:
		id_dict.setdefault(id_type, []).append(id_val)
	if len(mesh_ids) == 0:
	    continue
	id_to_mesh_ids[node] = mesh_ids
	for mesh_id in mesh_ids:
	    #mesh_id_to_type_to_ids[mesh_id] = id_dict
	    d = mesh_id_to_type_to_ids.setdefault(mesh_id, {})
	    for id_type, id_vals in id_dict.items():
		if id_type in d:
		    d[id_type] = list(set(id_vals + d[id_type]))
		else:
		    d[id_type] = id_vals
    return name_to_id, id_to_mesh_ids, mesh_id_to_type_to_ids


if __name__ == "__main__":
    main()

