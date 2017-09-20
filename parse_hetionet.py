import os, cPickle, parse_drugbank


def main():
    mapping_file = "/home/emre/data/indication/hetionet/hetnet/tsv/hetionet-v1.0-edges.sif"
    drug_to_do_ids = get_hetionet_mapping(mapping_file, metaedge="CtD")
    print drug_to_do_ids["DB00284"]
    drug_to_indications = get_hetionet_mapping(mapping_file)
    print drug_to_indications["DB00997"]
    return 


def get_drug_disease_mapping(drug_to_do_ids, do_to_mesh_ids, mesh_id_to_name, dump_file):
    """
    drug_to_do_ids from HetioNet
    do_to_mesh_ids from DO
    mesh_id_to_name mapping from UMLS (MSH | MH)
    """
    if dump_file is not None and os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_diseases = {} # (mesh_term, mesh_id, association score) 
    for drugbank_id, do_ids in drug_to_do_ids.iteritems():
        for do_id in do_ids:
	    if do_id in do_to_mesh_ids:
		for dui in do_to_mesh_ids[do_id]:
		    phenotype = mesh_id_to_name[dui].lower()
		    drug_to_diseases.setdefault(drugbank_id, set()).add((phenotype, dui, 1))
    if dump_file is not None:
	cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases


def get_hetionet_mapping(mapping_file, metaedge="CtD"):
    """ 
    Compound::DB00997	CtD Disease::DOID:184
    """
    drug_to_do_ids = {}
    f = open(mapping_file)
    header = f.readline()
    #m = len("Compound::")
    #n = len("Disease::")
    for line in f:
        words = line.strip().split("\t")
	edge_type = words[1]
	if edge_type != metaedge:
	    continue
	drug = words[0].split("::")[-1]
	do_id = words[2].split("::")[-1]
        drug_to_do_ids.setdefault(drug, set()).add(do_id)
    f.close()
    return drug_to_do_ids
    
	
if __name__ == "__main__":
    main()

