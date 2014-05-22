import urllib2, os, cPickle
from bs4 import BeautifulSoup
from diseasome.src import diseasome
from xml.etree.ElementTree import iterparse

def main():
    #drug = "metformin" 
    drug = "vitamin e" 
    #drug = "tesamorelin" #"egrifta" 
    print get_drug_treatment(drug)
    return


def get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id):
    disease_to_drugs = {}
    mesh_ids = set(phenotype_to_mesh_id.values())
    for drugbank_id, diseases in drug_to_diseases.iteritems():
	for disease, dui, val in diseases:
	    if dui in mesh_ids: # In the disease data set
		disease = disease.lower()
		#if phenotype_to_mesh_id[disease] != dui:
		#    print "Warning: Inconsistent dui", disease, phenotype_to_mesh_id[disease], dui
		disease_to_drugs.setdefault(disease, set()).add(drugbank_id)
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


def get_data(command, parameter):
    if command == "list":
        parameter = parameter.replace(" ", "+")
        txt = "search?conceptName=%s&kindName=DRUG_KIND" % (parameter)
    elif command == "get":
        txt = "allInfo/%s" % parameter 
    else:
        raise ValueError("Unknown command: " + command)
    url = 'http://rxnav.nlm.nih.gov/REST/Ndfrt/%s' % txt 
    #print url
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    return response


def get_drug_treatment(drug):
    #try:
    cid = get_drug_concept_id(drug)
    #except urllib2.HTTPError:
    if cid is None:
        #print "No info for", drug
    	return None
    #print cid
    diseases = get_concept_treatment_info(cid)
    diseases_mod = []
    for disease, cid, val in diseases:
        name, dui, cui = get_concept_mesh_info(cid)
	if dui is None: # Not a disease term, e.g., Menopause
	    continue
        diseases_mod.append((name, dui, val))
    return diseases_mod


def get_concept_mesh_info(cid):
    response = get_data("get", cid)
    soup = BeautifulSoup(response, "xml")
    name, dui, cui = None, None, None
    for tag in soup.find_all('propertyName'): # property
        #for tag_inner in tag.children:
        #    if tag_inner.name == "propertyName":
        #        word = tag_inner.string
        word = tag.string
	#print word
        if word == "MeSH_DUI":
            dui = str(tag.next_sibling.string) # tag_inner
        elif word == "MeSH_Name":
            name = tag.next_sibling.string.encode('ascii','ignore')
        elif word == "UMLS_CUI":
            cui = str(tag.next_sibling.string)
    return name, dui, cui


def get_concept_treatment_info(cid):
    response = get_data("get", cid)
    diseases = []
    # Beatiful Soup version
    soup = BeautifulSoup(response, "xml")
    for tag in soup.find_all('role'): #, recursive=False):
        for tag_inner in tag.children:
            #print tag_inner.name
            if tag_inner.name == "roleName":
                word = tag_inner.string
		if word.startswith("may_treat"):
                    flag = 0.5 
		elif word.startswith("treats"):
                    flag = 1 
                else:
                    flag = 0 
            elif tag_inner.name == "concept":
                for child in tag_inner.children:
                    #print child.name
                    if child.name == "conceptName":
                        disease = child.string
                    if child.name == "conceptKind":
                        if child.string.strip() != "DISEASE_KIND":
                            flag = 0
                    if child.name == "conceptNui":
                        disease_cid = child.string
                if flag > 0:
                    diseases.append((disease, disease_cid, flag))
    return diseases
    # ElemTree version
    context = iterparse(response, ["start", "end"])
    context = iter(context)
    event, root = context.next()
    NS = ""
    state_stack = [ root.tag ]
    flag = False
    for (event, elem) in context:
        if event == "start":
            state_stack.append(elem.tag)
            if elem.tag == NS+"role":
                flag = False
        elif event == "end":
            if elem.tag == NS+"roleName":
                if state_stack[-2] == NS+"role":
                    word = elem.text
                    if word.startswith("may_treat") or word.startswith("treats"):
                        flag = True
                    else:
                        flag = False
            if elem.tag == NS+"conceptName":
                if state_stack[-2] == NS+"concept" and state_stack[-3] == NS+"role":
                    if flag:
                        diseases.append(elem.text)
            elem.clear()
            state_stack.pop()
    root.clear()
    return diseases


def get_drug_concept_id(concept_name):
    cid = None
    response = get_data("list", concept_name)
    soup = BeautifulSoup(response, "xml")
    for tag in soup.find_all('concept'):
        for tag_inner in tag.children:
            #print tag_inner.name, tag_inner.contents
            if tag_inner.name == "conceptNui":
                values = tag_inner.contents
                if len(values) > 1:
                    print "Warning: concept has multiple ids", values
                cid = values[0]
    return cid


if __name__ == "__main__":
    main()

