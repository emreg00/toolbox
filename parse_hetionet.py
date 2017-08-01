import os, cPickle, parse_drugbank


def main():
    mapping_file = "/home/emre/arastirma/data/indication/hetionet/hetnet/tsv/hetionet-v1.0-edges.sif"
    drug_to_indications = get_hetionet_mapping(mapping_file)
    print drug_to_indications["DB00997"]
    return 


#!
def get_drug_disease_mapping(drug_to_do_ids, do_id_to_name, do_id_to_meshid, concept_id_to_mesh_id, mesh_id_to_name, dump_file):
    """
    name_to_cui_and_confidences mapping from MEDI
    name_to_drug & synonym mapping from DrugBank
    concept & mesh_id mapping from UMLS (MSH | MH)
    """
    if dump_file is not None and os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_diseases = {} # (mesh_id, mesh_term, association score) 
    for name, values in name_to_cui_and_confidences.iteritems():
        # Get drugbank id from name in the label
	drugbank_id, drugbank_name = parse_drugbank.get_drugbank_id_from_name(name, name_to_drug, synonym_to_drug)
        if drugbank_id is None:
            continue
	print "%s\t%s\t%s" % (drugbank_name, drugbank_id, name)
        for cui, val in values:
            if cui in concept_id_to_mesh_id:
                dui = concept_id_to_mesh_id[cui]
                phenotype = mesh_id_to_name[dui]
                drug_to_diseases.setdefault(drugbank_id, set()).add((phenotype, dui, val))
    if dump_file is not None:
	cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases

#!
def get_hetionet_mapping(mapping_file):
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
	if edge_type != "CtD":
	    continue
	drug = words[0].split("::")[-1]
	do_id = words[2].split("::")[-1]
        drug_to_do_ids.setdefault(drug, set()).add(do_id)
    f.close()
    return drug_to_do_ids
    
	
if __name__ == "__main__":
    main()

