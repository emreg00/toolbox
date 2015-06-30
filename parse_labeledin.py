import os, cPickle, re 


def main():
    spl_id = "9c983da1-a8e0-479e-8816-5d5668ea242e" 
    labeledin_file = "/home/emre/arastirma/data/drug/labeledin/LabeledIn_Structured_Results.txt"
    spl_id_to_cuis = get_labeledin_mapping(labeledin_file)
    print spl_id_to_cuis[spl_id]
    return 


def get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id):
    disease_to_drugs = {}
    mesh_id_to_phenotype = {}
    for phenotype, mesh_id in phenotype_to_mesh_id.items():
        mesh_id_to_phenotype[mesh_id] = phenotype
    for drugbank_id, diseases in drug_to_diseases.iteritems():
        for phenotype, dui, val in diseases:
            if val > 0:
                if dui in mesh_id_to_phenotype: # In the disease data set
                    disease = mesh_id_to_phenotype[dui].lower()
                    disease_to_drugs.setdefault(disease, set()).add(drugbank_id)
    return disease_to_drugs


def get_drug_disease_mapping(spl_id_to_cuis, spl_id_to_names, name_to_drug, synonym_to_drug, concept_id_to_mesh_id, mesh_id_to_name, dump_file):
    """
    spl_id_to_cuis mapping from LabeledIn
    spl_id_to_names mapping from DailyMed rxnorm
    name_to_drug & synonym mapping from DrugBank
    concept & mesh_id mapping from UMLS (MSH | MH)
    """
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_diseases = {} # (mesh_id, mesh_term, association score) 
    val = 1
    for spl_id, cuis in spl_id_to_cuis.iteritems():
        if spl_id not in spl_id_to_names:
            continue
	for name in spl_id_to_names[spl_id]:
	    name = name.lower()
	    # Get drugbank id from name in the label
	    drugbank_id = None
	    drugbank_name = None
	    for db_name, db_id in name_to_drug.iteritems():
		db_name = db_name.lower()
		exp = re.compile(r"\b%s\b" % db_name)
		m = exp.search(name)
		if m is None: 
		    continue
		#if drugbank_id is not None:
		#    print drugbank_id, db_id, name
		drugbank_id = db_id
		drugbank_name = db_name
	    if drugbank_id is None:
		for db_name, db_id in synonym_to_drug.iteritems():
		    db_name = db_name.lower()
                    try:
                        exp = re.compile(r"\b%s\b" % db_name)
                    except:
                        continue
                    m = exp.search(name)
		    if m is None: 
			continue
		    #if drugbank_id is not None:
		    #	print drugbank_id, db_id, name
		    drugbank_id = db_id
		    drugbank_name = db_name
	    if drugbank_id is not None:
		break
        if drugbank_id is None:
            continue
	print "%s\t%s\t%s\t%s" % (spl_id, drugbank_id, drugbank_name, name)
        for cui in cuis:
            if cui in concept_id_to_mesh_id:
                dui = concept_id_to_mesh_id[cui]
                phenotype = mesh_id_to_name[dui]
                drug_to_diseases.setdefault(drugbank_id, set()).add((phenotype, dui, val))
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases


def get_labeledin_mapping(mapping_file):
    """ 
    ID|SPL|CUI|RXCUIs...
    1749|2B1B7B5F-2F20-418C-B1AC-794C2EF1CE5E|C0002395;|135447|483071;|997216;|483071;|997216;
    """
    spl_id_to_cuis = {}
    f = open(mapping_file)
    for line in f:
        words = line.strip().split("|")
        spl_id = words[1].lower()
        cuis = [ word.strip() for word in words[2].split(";") if word.strip() != "" ]
        #rxcui_no_dosage = words[4].split(";")
        #rxcui = words[5].split(";")
        #rxcui_other_no_dosage = words[6].split(";")
        #rxcui_other = words[7].split(";")
        spl_id_to_cuis[spl_id] = cuis
    f.close()
    return spl_id_to_cuis
    
	
if __name__ == "__main__":
    main()

