import OBO

def main():
    medic_file = "../../data/ontology/CTD_diseases.obo"
    name_to_id, id_to_mesh_ids = get_medic_mesh_id_mapping(medic_file)
    return


def get_medic_mesh_id_mapping(medic_file):
    medic = OBO.OBO(medic_file, save_synonyms = True)
    name_to_id = {}
    id_to_mesh_ids = {}
    for node, data in medic.g.nodes(data=True):
	name = data['n']
	name_to_id[name] = node
	if 's' in data:
	    for name in data['s']:
		name_to_id[name] = node
	if node.startswith("MESH:D"):
	    id_to_mesh_ids[node] = [ node[5:] ]
	else:
	    for node2, type2 in medic.get_term_relations(node):
		if type2 == "is_a":
		    if node2.startswith("MESH:D"):
			id_to_mesh_ids.setdefault(node, []).append(node2[5:])
    return name_to_id, id_to_mesh_ids 


if __name__ == "__main__":
    main()

