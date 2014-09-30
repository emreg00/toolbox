import urllib2, os, cPickle, json

API_USER_KEY = None # change this value with custom API key

def main():
    drug = "metformin" 
    disease = "type 2 diabetes"
    #drug = "vitamin e" 
    #drug = "tesamorelin" #"egrifta" 
    values = get_drug_treatment(drug, disease)
    values.sort()
    print values
    return


def get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id):
    disease_to_drugs = {}
    #disease_id_to_drugs = {}
    mesh_ids = set(phenotype_to_mesh_id.values())
    for drugbank_id, diseases in drug_to_diseases.iteritems():
	for disease, dui, val in diseases:
	    if dui in mesh_ids: # In the disease data set
		disease = disease.lower()
		#if phenotype_to_mesh_id[disease] != dui:
		#    print "Warning: Inconsistent dui", disease, phenotype_to_mesh_id[disease], dui
		disease_to_drugs.setdefault(disease, set()).add(drugbank_id)
		#disease_id_to_drugs(dui, set()).add(drugbank_id)
    return disease_to_drugs


def get_drug_disease_mapping(selected_drugs, drug_to_name, drug_to_synonyms, dump_file):
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_diseases = {} # (mesh_id, mesh_term, may_treat) 
    #flag = False
    for drugbank_id in selected_drugs:
	#if drugbank_id == "DB04575":
	#    flag = True
	#if flag == False: 
	#    continue
	drug = drug_to_name[drugbank_id].lower()
	print drugbank_id, drug,
	diseases = get_drug_treatment(drug)
	if diseases is None:
	    if drugbank_id not in drug_to_synonyms:
		print
		continue
	    for synonym in drug_to_synonyms[drugbank_id]:
		drug = synonym.lower()
		diseases = get_drug_treatment(drug)
		if diseases is not None:
		    break
	if diseases is None or len(diseases) == 0:
	    print
	    continue
	print diseases
	for disease, dui, val in diseases:
	    disease = disease.lower()
	    drug_to_diseases.setdefault(drugbank_id, []).append((disease, dui, val))
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases


def get_data(command, parameter, parameter2=None):
    field_drug = "patient.drug.medicinalproduct"
    field_disease = "patient.drug.drugindication"
    field_effect = "patient.reaction.reactionmeddrapt"
    if command == "disease-drug":
        parameter = parameter.replace(" ", "+")
        txt = '%s:"%s"&count=%s.exact' % (field_disease, parameter, field_drug)
    elif command == "drug-disease":
        parameter = parameter.replace(" ", "+")
        txt = '%s:"%s"&count=%s.exact' % (field_drug, parameter, field_disease)
    elif command == "drug-disease-effect":
        parameter = parameter.replace(" ", "+")
        parameter2 = parameter2.replace(" ", "+")
        txt = '%s:"%s"+AND+%s:"%s"&count=%s.exact' % (field_drug, parameter, field_disease, parameter2, field_effect)
    else:
        raise ValueError("Unknown command: " + command)
    if API_USER_KEY is None:
        url = 'https://api.fda.gov/drug/event.json?search=%s&limit=1000' % txt 
    else:
        url = 'https://api.fda.gov/drug/event.json?api_key=%s&search=%s&limit=1000' % (API_USER_KEY, txt)
    print url
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    while True:
	try:
	    response = json.load(response)
	    break
	except:
	    print "Problem with response:", parameter, parameter2
	    response = urllib2.urlopen(req)
    return response


def get_drug_treatment(drug, disease):
    #try:
    response = get_data("drug-disease-effect", drug, disease)
    #except urllib2.HTTPError:
    if len(response) == 0:
        print "No info for", drug, disease
    	return []
    values = []
    for row in response["results"]:
        #indication = row["patient"]["drug"]["drugindication"]
        #effect = row["patient"]["reaction"]["reactionmeddrapt"]
        term = row["term"]
        count = row["count"]
        print term, count #indication, effect
        values.append((count, term))
    return values


if __name__ == "__main__":
    main()

