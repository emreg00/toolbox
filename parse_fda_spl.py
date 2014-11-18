import urllib2, os, cPickle, re, time
from bs4 import BeautifulSoup


OFFSET = 66000
LIMIT = 105000 #105000 # 40-50K seem depleted (based on the first 1000)

def main():
    #drug = "montelukast" 
    drug = "methotrexate"
    disease = "type 2 diabetes mellitus"
    #disease = "asthma"
    output_dir = "/home/emre/arastirma/data/drug/fda/spl/"
    #fetch_spl_data(output_dir, OFFSET, LIMIT)
    for spl in range(1000): #[ 75, 82, 392 ]: #[ 2, 3, 4, 5, 8, 4013, 64602, 92978 ]:
        try:
            name, indication = read_spl_data(output_dir + "%d.html" % spl)
        except:
            continue
        #print name, indication
    return


def get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id):
    disease_to_drugs = {}
    mesh_id_to_phenotype = {}
    for phenotype, mesh_id in phenotype_to_mesh_id.items():
        mesh_id_to_phenotype[mesh_id] = phenotype
    for drugbank_id, diseases in drug_to_diseases.iteritems():
        for phenotype, dui in diseases:
	    if dui in mesh_id_to_phenotype: # In the disease data set
		disease = mesh_id_to_phenotype[dui].lower()
                disease_to_drugs.setdefault(disease, set()).add(drugbank_id)
    return disease_to_drugs


def get_drug_disease_mapping(selected_drugs, drug_to_name, drug_to_synonyms, mesh_id_to_name, mesh_id_to_name_with_synonyms, dump_file):
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_diseases = {} # (mesh_id, mesh_term, n, n_max, ri) 
    exp = re.compile("-\d-")
    mesh_name_to_ids = {}
    for mesh_id, names in mesh_id_to_name_with_synonyms.iteritems():
        if mesh_id.startswith("Q"):
            continue
        for name in names:
            #name = name.decode('utf-8','ignore')
            name = name.replace(",", "").lower()
            mesh_name_to_ids.setdefault(name, []).append(mesh_id)
    not_found_in_mesh = set()
    modified_in_mesh = set()
    multiple_mesh_id = {}
    flag = False 
    for drugbank_id in selected_drugs:
	if drugbank_id == "DB00229":
	    flag = True
	if flag == False: 
	    continue
        # Find the most common drug name in FDA for the DrugBank drug
	names = [ drug_to_name[drugbank_id] ]
        if drugbank_id in drug_to_synonyms:
            for synonym in drug_to_synonyms[drugbank_id]:
                if synonym.find("[") != -1 or synonym.find("{") != -1:
                    continue
                m = exp.search(synonym)
                if m:
                    continue
                names.append(synonym)
	drug, n = choose_fda_drug_name(names)
        if drug is None: # No match in FDA api
            continue
	#print drugbank_id, drug, n
	# Get diseases for that drug in FDA
        diseases = get_diseases_for_drug(drug)
        #time.sleep(0.25) # 240 request / min limit
        if len(diseases) == 0:
            continue
        f = open(dump_file + ".txt", 'a')
        for n, disease in diseases:
            if n < N_MIN: 
                continue
            disease = disease.lower()
            phenotype = convert_fda_name_to_mesh(disease, mesh_name_to_ids)
            if phenotype is None:
                not_found_in_mesh.add(disease)
                continue
            duis = mesh_name_to_ids[phenotype]
            if len(duis) > 1:
                multiple_mesh_id[phenotype] = duis
                continue
            dui = duis[0]
            if dui not in mesh_id_to_name:
                continue
            phenotype_mod = mesh_id_to_name[dui]
            if phenotype != phenotype_mod: # matched to synonym or removed 's
                modified_in_mesh.add((phenotype, phenotype_mod))
            # API can not deal with x^s y / x's y such as crohn's disease and non-hodgkin's lymphoma
            idx = disease.find("^s")
            if idx == -1:
                idx = disease.find("'s")
            if idx != -1: 
                phenotype = disease[:idx]
            # Get efficacy for each disease
            values, values_eff = get_drug_treatment(drug, phenotype)
            time.sleep(0.3) # 240 request / min limit
            if values is None or len(values) == 0:
                continue
            z_max, count_max, term = values[0]
            z_ineff, count_ineff, z_adverse, count_adverse = values_eff
            # Safety reports provide multiple drug-disease pairs count_max is more reliable than n
            #print disease, phenotype, n, count_max, count_ineff, count_adverse
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (drugbank_id, drug_to_name[drugbank_id], drug, disease, phenotype_mod, dui, n, count_max, count_ineff, count_adverse))
            #drug_to_diseases.setdefault(drugbank_id, []).append((phenotype, dui, n, count_max, count_ineff, count_adverse))
        f.close()
    print "Not found in MeSH:", not_found_in_mesh
    print "Modified in MeSH:", modified_in_mesh
    print "Multiple id in MeSH:", multiple_mesh_id
    for line in open(dump_file + ".txt"):
        (drugbank_id, name, fda_name, fda_disease, phenotype, dui, n, count_max, count_ineff, count_adverse) = line.strip().split("\t")
        for mesh_id in dui.split(","):
            drug_to_diseases.setdefault(drugbank_id, []).append((phenotype, mesh_id, int(n), int(count_max), int(count_ineff), int(count_adverse)))
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases



def get_data(command, parameter):
    response = None
    if command == "drug":
        txt = '%s' % (parameter)
    else:
        raise ValueError("Unknown command: " + command)
    url = 'https://rm2.scinet.fda.gov/druglabel/rs/spl/by-id/%s/%s.html' % (txt, txt)
    #print url
    req = urllib2.Request(url)
    try:
        response = urllib2.urlopen(req)
    except:
        print "Problem with response:", parameter
    return response


def fetch_spl_data(output_dir, offset, limit):
    while offset < limit:
        offset += 1
        out_file = output_dir + "%s.html" % offset
        if not os.path.exists(out_file):
            result = get_data("drug", offset)
            if result is None:
                continue
            f = open(out_file, 'w')
            for row in result:
                f.write(row)
            f.close()
        #name, indication = read_spl_data(out_file)
        #print name, indication
    return #result


def read_spl_data(file_name):
    name = None
    indication = [] 
    html_doc = open(file_name)
    soup = BeautifulSoup(html_doc, "xml")
    for tag in soup.find_all('p', class_="DocumentTitle"):
        #print tag.name
        if name is not None:
            print "Multiple name:", name
        name = tag.strong.string.encode("ascii", "ignore")
        words = name.split("-")
        name = words[0].rstrip().lower()
    for tag in soup.find_all('h1'):
        #print tag.name 
        if tag.string is None:
            continue
        header = tag.string.encode().lower()
        #print header
        if header.find("indication") != -1 or header in ("uses", "use", "usage"): # "indications", "indications and usage"
            if header.find("contraindication") != -1:
                continue
            for tag_p in tag.find_all_next():
                try:
                    #print "II:", tag_p.name
                    if tag_p.name in "h1": #("h1", "h2"):
                        break
                except:
                    continue
                if tag_p.name == "p" or tag_p.name == "li":
                    txt = tag_p.get_text().replace("  ", "").replace("\n","").encode("ascii", "ignore")
                    if txt == "":
                        continue
                    indication.append(txt)
                    #if tag_p.string is None:
                    #    continue
            if len(indication) == 0:
                print "No indication:", name
                for tag_p in tag.next_siblings:
                    try:
                        #print "I:", tag_p.name
                        if tag_p.name in ("h1", "h2"):
                            break
                    except:
                        continue
                    #if tag_p.string is None:
                    #    continue
                    #txt = tag_p.string.replace("  ", "").replace("\n","")
                    txt = tag_p.get_text() 
                    #tag2 = tag.find_next("h1")
                    #txt2 = tag2.get_text()
                    #idx = txt.find(txt2)
                    #txt = txt[:idx]
                    txt = txt.replace("  ", "").replace("\n","")
                    txt = txt.encode("ascii", "ignore")
                    if txt == "":
                        continue
                    indication.append(txt)
            break
    return name, indication

    
	
if __name__ == "__main__":
    main()

