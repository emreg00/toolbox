import urllib2, os, cPickle
from bs4 import BeautifulSoup
import time

def main():
    #drug = "metformin" 
    #drug = "Egrifta" #"Tesamorelin"
    drug = "vitamin e" 
    print get_compound_to_mesh_concepts(drug, 0.05)
    return


def get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id, cutoff_qvalue = 0.05):
    disease_to_drugs = {}
    mesh_ids = set(phenotype_to_mesh_id.values())
    for drugbank_id, diseases in drug_to_diseases.iteritems():
	for disease, dui, val in diseases:
	    if dui in mesh_ids: # In the disease data set
		disease = disease.lower()
		#if phenotype_to_mesh_id[disease] != dui:
		#    print "Warning: Inconsistent dui", disease, phenotype_to_mesh_id[disease], dui
		if float(val) <= cutoff_qvalue:
		    disease_to_drugs.setdefault(disease, set()).add(drugbank_id)
    return disease_to_drugs


def get_drug_disease_mapping(selected_drugs, drug_to_name, drug_to_synonyms, mesh_id_to_name, dump_file):
    if os.path.exists(dump_file):
    	drug_to_diseases = cPickle.load(open(dump_file))
    	return drug_to_diseases 
    mesh_ids = set(mesh_id_to_name.keys())
    drug_to_diseases = {}
    #flag = False
    for drugbank_id in selected_drugs:
	#if drugbank_id == "DB00956":
	#    flag = True
	#if flag == False:
	#    continue
	drug = drug_to_name[drugbank_id].lower()
	print drugbank_id, drug,
	diseases = get_compound_to_mesh_concepts(drug) #, 1e-10) 
	if len(diseases) == 0:
	    if drugbank_id not in drug_to_synonyms:
		print drugbank_id
		continue
	    for synonym in drug_to_synonyms[drugbank_id]:
		drug = synonym.lower()
		diseases = get_compound_to_mesh_concepts(drug) #, 1e-10) 
		if len(diseases) > 0:
		    break
	if len(diseases) == 0:
	    print
	    continue
        f = open(dump_file + ".txt", 'a')
	diseases_mod = []
	for phenotype, dui, q_value in diseases:
	    if dui in mesh_ids: 
		#drug_to_diseases.setdefault(drugbank_id, []).append((phenotype, dui, q_value))
		diseases_mod.append((phenotype, dui, q_value))
		f.write("%s\t%s\t%s\t%e\n" % (drugbank_id, phenotype, dui, q_value))
	    #else:
	    #	print "\nNot in UMLS:", dui, phenotype
        f.close()
	print diseases_mod
    for line in open(dump_file + ".txt"):
        (drugbank_id, phenotype, dui, q_value) = line.strip().split("\t")
        drug_to_diseases.setdefault(drugbank_id, []).append((phenotype, dui, q_value))
    cPickle.dump(drug_to_diseases, open(dump_file,'w'))
    return drug_to_diseases


def get_data(command, parameter):
    """
    command: compound | mesh
    """
    parameter = parameter.replace(" ", "+")
    url = 'http://metab2mesh.ncibi.org/fetch?%s=%s&limit=9999&publimit=1' % (command, parameter)
    #print url
    req = urllib2.Request(url)
    req.add_header('User-agent', 'Mozilla/5.0 (Linux i686)')
    while True:
	try:
	    response = urllib2.urlopen(req)
	except urllib2.URLError:
            time.sleep(10)
	    continue
	break
    return response


def get_compound_to_mesh_concepts(drug, cutoff_qvalue = 0.05):
    response = get_data("compound", drug)
    phenotypes = []
    # Beatiful Soup version
    soup = BeautifulSoup(response, "xml")
    for tag in soup.find_all('Result'): #, recursive=False):
        phenotype = None
        for tag_inner in tag.children:
            #print tag_inner.name
            if tag_inner.name == "MeSH":
                for child in tag_inner.Descriptor.children:
                    #print child.name
                    if child.name == "Name":
                        phenotype = child.string.encode("ascii", "ignore")
                    if child.name == "Identifier":
                        phenotype_cid = str(child.string)
            elif tag_inner.name == "FisherExact":
                p_value = float(tag_inner.string)
            elif tag_inner.name == "Q-Value":
                q_value = float(tag_inner.string)
        if phenotype is not None:
            if q_value <= cutoff_qvalue:
                phenotypes.append((phenotype, phenotype_cid, q_value))
    return phenotypes


if __name__ == "__main__":
    main()

