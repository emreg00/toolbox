import os, cPickle, parse_drugbank


def main():
    name = "Mesna"
    mapping_file = "/home/emre/arastirma/data/drug/medi/MEDI_01212013_UMLS.csv"
    name_to_cui_and_confidences = get_medi_mapping(mapping_file)
    print name_to_cui_and_confidences[name]
    return 


def get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id):
    disease_to_drugs = {}
    mesh_id_to_phenotype = {}
    for phenotype, mesh_id in phenotype_to_mesh_id.items():
        mesh_id_to_phenotype[mesh_id] = phenotype
    for drugbank_id, diseases in drug_to_diseases.iteritems():
        for phenotype, dui, val in diseases:
            if val > 0.5:
                if dui in mesh_id_to_phenotype: # In the disease data set
                    disease = mesh_id_to_phenotype[dui].lower()
                    disease_to_drugs.setdefault(disease, set()).add(drugbank_id)
    return disease_to_drugs


def get_drug_disease_mapping(name_to_cui_and_confidences, name_to_drug, synonym_to_drug, concept_id_to_mesh_id, mesh_id_to_name, dump_file):
    """
    name_to_cui_and_confidences mapping from MEDI
    name_to_drug & synonym mapping from DrugBank
    concept & mesh_id mapping from UMLS (MSH | MH)
    """
    if os.path.exists(dump_file):
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
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases


def get_medi_mapping(mapping_file):
    """ 
    RXCUI_IN,DRUG_DESC,UMLS_CUI,ICD9,INDICATION_DESCRIPTION,MENTIONEDBYRESOURCES,HIGHPRECISIONSUBSET,POSSIBLE_LABEL_USE
    44,Mesna,C0018965,599.7,Hematuria,1,0,0
    """
    name_to_cui_and_confidences = {}
    f = open(mapping_file)
    header = f.readline()
    for line in f:
        words = line.strip().split(",")
        name = words[1]
        cui = words[2]
        in_hps = words[6]
        confidence = 0.5
        if in_hps == "1":
            confidence = 1 
        name_to_cui_and_confidences.setdefault(name, []).append((cui, confidence))
    f.close()
    return name_to_cui_and_confidences
    
	
if __name__ == "__main__":
    main()

