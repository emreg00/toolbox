import urllib2, os, cPickle, re, time
from bs4 import BeautifulSoup
from toolbox import text_utilities


OFFSET = 1
LIMIT = 111000 # 40K is depleted 

def main():
    #drug = "montelukast" 
    #drug = "methotrexate"
    #disease = "type 2 diabetes mellitus"
    #disease = "asthma"
    output_dir = "/home/emre/arastirma/data/drug/fda/spl/"
    for spl in (1083, 1156, 3236, 3623, 12235):
        name, indication, contraindication, warning = read_spl_data(output_dir + "%d.html" % spl)
        print spl, name
        print indication
        print contraindication
        print warning
    return 
    #fetch_spl_data(output_dir, OFFSET, LIMIT)
    for spl in range(19000,20000): #[ 2, 3, 4, 5, 8, 4013, 64602, 92978 ]: #[ 75, 82, 392, 10014 ]: #
        print spl
        try:
            name, indication, contraindication, warning = read_spl_data(output_dir + "%d.html" % spl)
        except:
            continue
        print name #, indication
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


def get_drug_disease_mapping(output_dir, selected_drugs, name_to_drug, synonym_to_drug, mesh_id_to_name, mesh_id_to_name_with_synonyms, negex_file, dump_file):
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_diseases = {} # (mesh_id, mesh_term, non-symptomaticy score) 
    mesh_name_to_id = {}
    abbreviation_to_id = {}
    # Mesh mapping already filtered for those that are disease terms
    for mesh_id, names in mesh_id_to_name_with_synonyms.iteritems():
        for name in names:
            # Take into account abbreviations (<= 4 letter) and the case of AIDS
            if name.isupper() and len(name) > 1:
		abbreviation_to_id[name] = mesh_id
            else:
                name = name.lower()
		# Remove the final s (plural) and match later with an additional s
		if name.endswith('s'):
		    name = name[:-1] 
		#name = " " + name 
		#name = name.decode('utf-8','ignore')
		for name_mod in [ name, name.replace(",", ""), name.replace("-", " "), name.replace(",", "").replace("-", " ") ]:
		    mesh_name_to_id[name_mod] = mesh_id
    # Get keywords / negex for text matching
    negex_rules = text_utilities.get_negex_rules(negex_file)
    flag = False 
    for spl in xrange(LIMIT):
	if spl == OFFSET:
	    flag = True
	if flag == False: 
	    continue
        try:
            name, indication, contraindication, warning = read_spl_data(output_dir + "%d.html" % spl)
            #print "SPL:", spl #, name
            if len(warning) == 0: # Skip potentially problematic label
                continue
        except:
            continue
        # Get drugbank id from name in the label
        drugbank_id = None
        if name in name_to_drug:
            drugbank_id = name_to_drug[name]
        elif name in synonym_to_drug:
            drugbank_id =  synonym_to_drug[name]
        else:
            continue
        if selected_drugs is not None and drugbank_id not in selected_drugs: # Wont happen since name mapping used only selected_drugs
            #print "Not in selected:", drugbank_id
            continue
	print spl, drugbank_id, name
        # Sentencify
        for idx, values in enumerate([indication, contraindication, warning]):
            indications = []
            for txt in values:
                for sentence in txt.lower().split("."):
                    indications.append(sentence)
            if len(indications) == 0:
                continue
            # Match the indication to mesh keywords
            for mesh_name, dui in mesh_name_to_id.iteritems():
                exp = re.compile(r"\b%ss{,1}\b" % mesh_name)
                for sentence in indications:
                    # Look for mesh term 
                    m = exp.search(sentence)
                    if m is None: 
                        continue
                    # Was removing negative sentences from indication: not / except / no / inappropriate, now assign -1 score
                    val = get_value_of_association(mesh_name, sentence, negex_rules, idx)
                    #if dui not in mesh_id_to_name: 
                    #    continue
                    phenotype = mesh_id_to_name[dui]
                    #if val != 1:
                    #print "A/S/N:", mesh_name, val, phenotype, dui, sentence
                    drug_to_diseases.setdefault(drugbank_id, set()).add((phenotype, dui, val))
            # Handle the abbreviations separately (keep the indication "upper" and match against it)
            for mesh_name, dui in abbreviation_to_id.iteritems():
                exp = re.compile(r"\b%ss{,1}\b" % mesh_name)
                for sentence in values:
                    m = exp.search(sentence)
                    if m is None:
                        continue
                    val = get_value_of_association(mesh_name, sentence, negex_rules, idx)
                    phenotype = mesh_id_to_name[dui]
                    drug_to_diseases.setdefault(drugbank_id, set()).add((phenotype, dui, val))
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases


def get_value_of_association(mesh_name, sentence, negex_rules, idx):
    val = 1
    # Symptomatic cases: protect / maintain / manage(ment) / symptom / relie(f) - relie(ve) / palliati(ve) - palliati(on) / alleviate
    if idx == 0: # indication
        symptomatic, i = text_utilities.is_symptomatic(sentence)
        if symptomatic:
            if i == 0:
                val = 0.8 # protection
            elif i == 1:
                val = 0.7 # maintain / maintenance
            elif i == 2:
                val = 0.6 # manage(ment)
            else:
                val = 0.5
        negative = text_utilities.is_negated(sentence, mesh_name, negex_rules)
        negative2 = text_utilities.is_negated(sentence, mesh_name, None) 
        if negative and negative2:
            val = -0.7
        elif negative:
            val = -0.6
        elif negative2:
            val = -0.5
        #if negative != negative2: 
        #    print "N:", mesh_name, negative, negative2, sentence
    elif idx == 1: # contraindication
        val = -1
    elif idx == 2: # warning
        val = -0.9
    else:
        raise ValueError("Unexpected index!")
    return val


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
    contraindication = [] 
    warning = []
    html_doc = open(file_name)
    soup = BeautifulSoup(html_doc, "xml")
    for tag in soup.find_all('p', class_="DocumentTitle"):
        #print tag.name
        if name is not None:
            print "Multiple name:", name
        name = tag.strong.string.encode("ascii", "ignore")
        words = name.split(" - ")
        name = words[0].strip().lower()
    for tag in soup.find_all('h1'):
        #print tag.name 
        if tag.string is None:
            continue
        header = tag.string.encode().strip().lower()
        #print header
        flag = None
        if header.find("contraindication") != -1:
            flag = "contraindication"
        elif header.find("indication") != -1 or header in ("uses", "use", "usage"): # "indications", "indications and usage"
            flag = "indication"
        if header.find("warning") != -1 and header.find("boxed") == -1:
            flag = "warning"
        if flag is not None:
            for tag_p in tag.find_all_next():
                try:
                    #print "II:", tag_p.name
                    if tag_p.name == "h1": # in ("h1", "h2"):
                        break
                except:
                    continue
                if tag_p.name == "p" or tag_p.name == "li":
                    txt = tag_p.get_text().strip()
                    if txt == "":
                        continue
                    tag2 = tag.find_next("h1")
                    txt2 = tag2.get_text().strip()
                    idx = txt.find(txt2)
                    if idx != -1:
                        txt = txt[:idx]
                    txt = " ".join(txt.split()) 
                    txt = txt.encode("ascii", "ignore")
                    if txt == "":
                        continue
                    if flag == "indication":
                        indication.append(txt)
                    elif flag == "contraindication":
                        contraindication.append(txt)
                    elif flag == "warning":
                        warning.append(txt)
            if (flag == "indication" and len(indication) == 0) or (flag == "contraindication" and len(contraindication) == 0) or (flag == "warning" and len(warning) == 0):
                # The indications are not encapsulated in <p>
                #print "No indication:", name
                for i, tag_p in enumerate(tag.next_siblings):
                    try:
                        #print "I:", tag_p.name
                        #print tag_p.string
                        if i == 0 and tag_p.string is not None:
                            txt = tag_p.string
                            txt = " ".join(txt.split()) 
                            txt = txt.encode("ascii", "ignore")
                            if txt == "":
                                continue
                            if flag == "indication":
                                indication.append(txt)
                            elif flag == "contraindication":
                                contraindication.append(txt)
                            elif flag == "warning":
                                warning.append(txt)
                        if tag_p.name == "h1":  
                            break
                    except:
                        continue
                    #txt = tag_p.string.replace("  ", " ").replace("\t"," ").replace("\n","")
                    txt = tag_p.get_text().strip()
                    if txt == "":
                        continue
                    tag2 = tag.find_next("h1")
                    txt2 = tag2.get_text().strip()
                    idx = txt.find(txt2)
                    if idx != -1:
                        txt = txt[:idx]
                    txt = " ".join(txt.split()) 
                    txt = txt.encode("ascii", "ignore")
                    if txt == "":
                        continue
                    if flag == "indication":
                        indication.append(txt)
                    elif flag == "contraindication":
                        contraindication.append(txt)
                    elif flag == "warning":
                        warning.append(txt)
    return name, indication, contraindication, warning

    
	
if __name__ == "__main__":
    main()

