import urllib2, os, cPickle, json, re, time
import configuration, text_utilities, stat_utilities 

CONFIG = configuration.Configuration() 

try:
    API_USER_KEY = CONFIG.get("FDA_API_KEY") # change this value with custom API key
except:
    print "API key not found!"
    API_USER_KEY = None
LIMIT = 100
FIELD_DRUG = "patient.drug.medicinalproduct"
FIELD_DISEASE = "patient.drug.drugindication"
FIELD_EFFECT = "patient.reaction.reactionmeddrapt"

def main():
    #drug = "montelukast" 
    #drug = "vitamin e" # "vit. e"
    #drug = "donepezil"
    drug = "methotrexate"
    disease = "type 2 diabetes"
    #disease = "alzheimer"
    #disease = "asthma"
    #disease = "acute lymphocytic leukaemia"
    #disease = "rheumatoid arthritis"
    #disease = "bone sarcoma"
    condition = None #"drug ineffective"
    mesh_names = [ "Diabetes Mellitus, Type 2", "Osteopenia", "Bone Diseases, Metabolic", "Pulmonary Disease, Chronic Obstructive"]
    print convert_fda_name_to_mesh(disease, mesh_names)
    return
    print get_counts_for_drug(drug, disease, condition)
    #d = get_counts_from_data(drug, disease, condition)
    #print d, len(d)
    values, values_eff = get_drug_treatment(drug, disease)
    print map(lambda x: "%.2f(%d) %s" % x, values[:5])
    print values_eff
    #values = get_drugs_for_disease(disease)
    #print map(lambda x: "%.2f(%d) %s" % x, values[:20])
    #values = get_diseases_for_drug(drug)
    #print map(lambda x: "%.2f(%d) %s" % x, values[:20])
    return


def get_disease_specific_drugs(drug_to_diseases, phenotype_to_mesh_id):
    disease_to_drugs = {}
    mesh_id_to_phenotype = {}
    for phenotype, mesh_id in phenotype_to_mesh_id.items():
        mesh_id_to_phenotype[mesh_id] = phenotype
    for drugbank_id, diseases in drug_to_diseases.iteritems():
	for disease, dui, val in diseases:
	    if dui in mesh_id_to_phenotype: # In the disease data set
		disease = mesh_id_to_phenotype[dui].lower()
                disease_to_drugs.setdefault(disease, set()).add(drugbank_id)
    return disease_to_drugs


def get_drug_disease_mapping(selected_drugs, drug_to_name, drug_to_synonyms, mesh_id_to_name, dump_file):
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    drug_to_diseases = {} # (mesh_id, mesh_term, n, n_max, ri) 
    exp = re.compile("-\d-")
    #flag = False
    mesh_name_to_ids = {}
    for mesh_id, names in mesh_id_to_name.iteritems():
        for name in names:
            mesh_name_to_ids.setdefault(name, set()).add(mesh_id)
    for drugbank_id in selected_drugs:
	#if drugbank_id == "DB04575":
	#    flag = True
	#if flag == False: 
	#    continue
	names = [ drug_to_name[drugbank_id] ]
	for synonym in drug_to_synonyms[drugbank_id]:
	    if synonym.find("[") != -1 or synonym.find("{") != -1:
		continue
	    m = exp.search(synonym)
	    if m:
		continue
	    names.append(synonym)
	drug = choose_fda_drug_name(names)
	print drugbank_id, drug
	# Get diseases
	diseases = get_diseases_for_drug(drug)
        if len(diseases) == 0:
            continue
        f = open(dump_file + ".txt", 'a')
        not_found_in_mesh = set()
        for n, disease in diseases:
            phenotype = convert_fda_name_to_mesh(disease, mesh_name_to_ids.keys())
            if phenotype is None:
                not_found_in_mesh.add(disease)
                continue
            # Get efficacy for each disease
            values, values_eff = get_drug_treatment(drug, disease)
            if values is None or len(values) == 0:
                continue
            z_max, count_max, term = values[0]
            z_ineff, count_ineff, z_adverse, count_adverse = values_eff
            print disease, phenotype, n, count_max, count_ineff, count_adverse
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (drugbank_id, drug_to_name[drugbank_id], drug, disease, phenotype, ",".join(mesh_name_to_ids[phenotype])))
        f.close()
            #! for dui in mesh_name_to_ids[phenotype]:
            #!    drug_to_diseases.setdefault(drugbank_id, []).append((phenotype, dui, n, count_max, count_ineff, count_adverse))
    print "Not found in MeSH:", not_found_in_mesh
    for line in open(dump_file + ".txt"):
        (drugbank_id, phenotype, dui, n, count_max, count_ineff, count_adverse) = line.strip().split("\t")
        for mesh_id in dui.split(","):
            drug_to_diseases.setdefault(drugbank_id, []).append((phenotype, mesh_id, n, count_max, count_ineff, count_adverse))
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    return drug_to_diseases


def convert_mesh_name_to_fda_name(disease):
    disease_to_term = { "arrhythmias, cardiac": "arrhythmia", "colitis, ulcerative": "colitis ulcerative", "lung diseases, obstructive": "chronic obstructive pulmonary disease", "sarcoma": "bone sarcoma", "liver cirrhosis": "hepatic cirrhosis", "liver cirrhosis, biliary": "biliary cirrhosis primary", "anemia, hemolytic": "anemia", "coronary artery disease": "coronary artery disease", "varicose veins": "varicose vein", "blood coagulation disorders": "disseminated intravascular coagulation", "mycobacterium infections": "mycobacterium avium complex infection" } 
    if disease in disease_to_term:
        disease_mod = disease_to_term[disease] 
    # Chop "disease" (alzheimer, crohn),
    elif disease.endswith(" disease"):
        disease_mod = disease[:-len(" disease")]
    # Replace "neoplasms" with "cancer"
    elif disease.endswith(" neoplasms"):
        disease_mod = disease[:-len(" neoplasms")]
        disease_mod += " cancer"
    # Get rid of the "s" for disorders / diseases
    elif disease.endswith(" disorders"):
        disease_mod = disease[:-1]
    elif disease.endswith(" diseases"):
        disease_mod = disease[:-1]
    # Reverse the order for "," ("diabetes mellitus, type 2", "type 2 diabetes mellitus") 
    elif "," in disease:
        words = disease.split(", ")
        words.reverse()
        disease_mod = " ".join(words)
    else:
        disease_mod = disease.lower()
    return disease_mod


def convert_fda_name_to_mesh(disease, mesh_names):
    # Get words skipping disease / disorder / syndrome / plural / 's
    disease = disease.lower()
    values = text_utilities.tokenize_disease_name(disease, exact=False)
    #values_mod = []
    #for word in values:
        #if word == "neoplasm":
        #    word = "cancer"
        #elif word == "ii":
        #    word = "2"
        #elif word == "i":
        #    word = "1"
    #    values_mod.append(word)
    phenotype = None
    val_and_phenotypes = []
    for mesh_name in mesh_names:
        try: #!
            mesh_name = mesh_name.decode('utf-8','ignore')
        except:
            print mesh_name
        val = sum([ mesh_name.lower().find(word.strip()) != -1 for word in values_mod ])
        #print mesh_name, val
        if val > len(values_mod) / 2.0:
            #print mesh_name, disease
            #phenotype = mesh_name
            val_and_phenotypes.append((float(val)/len(mesh_name.split()), mesh_name))
    #print values, values_mod
    #print val_and_phenotypes
    if len(val_and_phenotypes) > 0:
        val_and_phenotypes.sort()
        phenotype = val_and_phenotypes[-1][1]
    return phenotype


def get_data_helper(command, parameter, parameter2=None, parameter_effect=None, skip=0):
    parameter = parameter.replace(" ", "+").replace("-", "+")
    if command == "drug":
        txt = '%s:"%s"' % (FIELD_DRUG, parameter)
    elif command == "disease":
        txt = '%s:"%s"' % (FIELD_DISEASE, parameter)
    elif command == "drug-disease":
	assert parameter2 is not None
	parameter2 = parameter2.replace(" ", "+").replace("-", "+")
	txt = '%s:"%s"+AND+%s:"%s"' % (FIELD_DRUG, parameter, FIELD_DISEASE, parameter2)
    elif command == "drug-disease-effect":
	assert (parameter2 is not None and parameter_effect is not None)
	parameter2 = parameter2.replace(" ", "+").replace("-", "+")
	parameter_effect = parameter_effect.replace(" ", "+").replace("-", "+")
	txt = '%s:"%s"+AND+%s:"%s"+AND+%s:"%s"' % (FIELD_DRUG, parameter, FIELD_DISEASE, parameter2, FIELD_EFFECT, parameter_effect)
    else:
        raise ValueError("Unknown command: " + command)
    if API_USER_KEY is None:
        url = 'https://api.fda.gov/drug/event.json?search=%s&limit=%d&skip=%d' % (txt, LIMIT, skip)
    else:
        url = 'https://api.fda.gov/drug/event.json?api_key=%s&search=%s&limit=%d&skip=%d' % (API_USER_KEY, txt, LIMIT, skip)
    #print url
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    while True:
	try:
	    response = json.load(response)
	    break
	except:
	    print "Problem with response:", parameter, parameter2
	    response = urllib2.urlopen(req)
    #n = int(response["meta"]["results"]["total"])
    return response


def get_data(command, parameter, parameter2=None, parameter_effect=None):
    offset = 0
    limit = LIMIT
    result = []
    while True:
	result2 = get_data_helper(command, parameter, parameter2, parameter_effect, skip=offset)
	print offset, len(result2["results"])
	result += result2["results"]
	offset += limit
	if len(result2["results"]) < limit:
	    break
    return result


def choose_fda_drug_name(names):
    values = []
    for name in names:
        name = name.lower()
        #if name.startswith("vitamin"):
        #    name = name.replace(" ", "+")
        #else:
        #    words = name.split(" ")
        #    if len(words) > 1: 
        #        #print "Chopping drug name", name
        #        name = words[0]
        #    if name in ("compound", "dr.", "salicylate", "sodium"):
        #        continue
        N, n, M, k = get_counts_for_drug(name, None, None)
        time.sleep(0.3) # 240 request / min limit
        values.append((N, name))
    values.sort()
    name = values[-1][1]
    return name


def get_counts_from_data(drug, disease, condition=None):
    if condition is None:
	condition_to_count = {}
	values = get_data("drug-disease", drug, disease)
	flag_drug = False
	flag_disease = False
	for row in values:
	    for row_inner in row["patient"]["reaction"]:
		if "reactionmeddrapt" not in row_inner:
		    continue
		condition = row_inner["reactionmeddrapt"].lower()
		i = condition_to_count.setdefault(condition, 0)
		condition_to_count[condition] = i + 1
	return condition_to_count
    values = get_data("drug-disease-effect", drug, disease, condition)
    i = 0
    flag_drug = False
    flag_disease = False
    for row in values:
	#print row["safetyreportid"]
	for row_inner in row["patient"]["drug"]:
	    if row_inner["medicinalproduct"].lower().find(drug) != -1:
		flag_drug = True
		#print row_inner["medicinalproduct"]
		if "drugindication" in row_inner:
		    if row_inner["drugindication"].lower() == disease:
			flag_disease = True
			#print row_inner["drugindication"]
	for row_inner in row["patient"]["reaction"]:
	    if "reactionmeddrapt" not in row_inner:
		continue
	    if row_inner["reactionmeddrapt"].lower() == condition:
		if flag_drug and flag_disease:
		    i += 1
		#print row_inner["reactionmeddrapt"]
    return i


def get_counts(command, parameter, parameter2=None, parameter_effect=None):
    parameter_org = parameter
    parameter = parameter.replace(" ", "+").replace("-", "+")
    if parameter2 is not None:
	parameter2_org = parameter2
        parameter2 = parameter2.replace(" ", "+").replace("-", "+")
    if parameter_effect is not None:
        parameter_effect = parameter_effect.replace(" ", "+").replace("-", "+")
    if command == "drug": # number of safety reports for that drug
        txt = '%s:"%s"&count=%s.exact' % (FIELD_DRUG, parameter, FIELD_DRUG)
    elif command == "disease": # number of safety reports for that disease
        txt = '%s:"%s"&count=%s.exact' % (FIELD_DISEASE, parameter, FIELD_DISEASE)
    elif command == "drug-disease": # number of safety reports for that drug and disease pair
	assert parameter2 is not None
	txt = '%s:"%s"+AND+%s:"%s"&count=%s.exact' % (FIELD_DRUG, parameter, FIELD_DISEASE, parameter2, FIELD_DRUG)
    elif command == "disease-drug":
	assert parameter2 is not None
        txt = '%s:"%s"+AND+%s:"%s"&count=%s.exact' % (FIELD_DISEASE, parameter2, FIELD_DRUG, parameter, FIELD_DISEASE)
    elif command == "drug-effect": # number of safety reports for that drug and reaction pair
	assert parameter_effect is not None
	txt = '%s:"%s"+AND+%s:"%s"&count=%s.exact' % (FIELD_DRUG, parameter, FIELD_EFFECT, parameter_effect, FIELD_DRUG)
    elif command == "disease-effect": # number of safety reports for that disease and reaction pair
	assert parameter_effect is not None
	txt = '%s:"%s"+AND+%s:"%s"&count=%s.exact' % (FIELD_DISEASE, parameter, FIELD_EFFECT, parameter_effect, FIELD_DISEASE)
    elif command == "drug-disease-effect" or command == "disease-drug-effect": # number of safety reports for that drug, disease and reaction triple
	assert (parameter2 is not None and parameter_effect is not None)
	if command == "drug-disease-effect":
	    txt = '%s:"%s"+AND+%s:"%s"+AND+%s:"%s"&count=%s.exact' % (FIELD_DRUG, parameter, FIELD_DISEASE, parameter2, FIELD_EFFECT, parameter_effect, FIELD_DRUG)
	elif command == "disease-drug-effect":
	    txt = '%s:"%s"+AND+%s:"%s"+AND+%s:"%s"&count=%s.exact' % (FIELD_DRUG, parameter, FIELD_DISEASE, parameter2, FIELD_EFFECT, parameter_effect, FIELD_DISEASE)
    elif command == "drug-disease2": # returns all diseases
        txt = '%s:"%s"&count=%s.exact' % (FIELD_DRUG, parameter, FIELD_DISEASE)
    elif command == "disease-drug2": # returns all drugs
        txt = '%s:"%s"&count=%s.exact' % (FIELD_DISEASE, parameter, FIELD_DRUG)
    elif command == "drug-disease-effect2": # returns all reactions and their counts
	assert parameter2 is not None
        txt = '%s:"%s"+AND+%s:"%s"&count=%s.exact' % (FIELD_DRUG, parameter, FIELD_DISEASE, parameter2, FIELD_EFFECT)
    else:
        raise ValueError("Unknown command: " + command)
    if API_USER_KEY is None:
        url = 'https://api.fda.gov/drug/event.json?search=%s&limit=%d' % (txt, 10*LIMIT)
    else:
        url = 'https://api.fda.gov/drug/event.json?api_key=%s&search=%s&limit=%d' % (API_USER_KEY, txt, 10*LIMIT)
    #print url
    n = None
    req = urllib2.Request(url)
    try:
	response = urllib2.urlopen(req)
    except urllib2.HTTPError:
	if parameter_effect is not None:
	    n = 0
	#print "No info for", parameter, parameter2, parameter_effect
	if parameter2 is not None or parameter_effect is not None:
	    print url
        return n
    while True:
	try:
	    response = json.load(response)
	    break
	except:
	    print "Problem with response:", parameter, parameter2
	    response = urllib2.urlopen(req)
    if command.endswith("2"):
	return response["results"]
    val = parameter_org.lower() #parameter.lower().replace("+", " ")
    if command in ("disease-drug", "disease-drug-effect"):
	val = parameter2_org.lower() #parameter2.lower().replace("+", " ")
    for row in response["results"]:
	# note that 's are ^s in the results
	if row["term"].lower().find(val) != -1:
	    #print row["term"].lower()
	    n = int(row["count"])
	    break
    return n


def z_scorize_counts(count_term_pairs):
    #values = []
    #for count, term in count_term_pairs:
    #	if count < 2: 
    #	    continue
    #	values.append((count, term))
    #count_term_pairs = values
    values = []
    m, s = stat_utilities.calc_mean_and_sigma(zip(*count_term_pairs)[0])
    for count, term in count_term_pairs:
        val = count - m
        if s != 0:
            val /= s
        values.append((val, count, term))
    values.sort()
    values.reverse()
    return values


def get_efficacy_values(values):
    z_ineff, count_ineff = 0, 0
    z_adverse, count_adverse = 0, 0
    for z, count, term in values:
        if term == "CONDITION AGGRAVATED":
            z_adverse += z
            count_adverse += count
        elif term == "DRUG INEFFECTIVE":
            z_ineff += z
            count_ineff += count
    return z_ineff, count_ineff, z_adverse, count_adverse


def get_counts_for_drug(drug, disease=None, condition=None):
    command = "drug"
    N = get_counts(command, drug)
    n, k, M = None, None, None
    if disease is not None:
	command = "drug-disease"
	n = get_counts(command, drug, disease)
	if condition is not None:
	    command = "drug-disease-effect"
	    k = get_counts(command, drug, disease, condition)
    if condition is not None:
	command = "drug-effect"
	M = get_counts(command, drug, None, condition)
    return N, n, M, k


def get_counts_for_disease(drug, disease=None, condition=None):
    command = "disease"
    N = get_counts(command, disease)
    n, k, M = None, None, None
    if disease is not None:
	command = "drug-disease" #"disease-drug"
	n = get_counts(command, drug, disease)
	if condition is not None:
	    command = "drug-disease-effect" #"disease-drug-effect"
	    k = get_counts(command, drug, disease, condition)
    if condition is not None:
	command = "disease-effect"
	M = get_counts(command, disease, None, condition)
    return N, n, M, k


def get_counts_for_drug_and_disease(drug, disease, condition=None):
    n, k = None, None
    command = "disease-drug"
    n = get_counts(command, drug, disease)
    if condition is not None:
	command = "disease-drug-effect"
	k = get_counts(command, drug, disease, condition)
    return n, k


def get_drug_treatment(drug, disease):
    response = get_counts("drug-disease-effect2", drug, disease)
    if response is None:
	return None, None
    values = []
    for row in response:
        #indication = row["patient"]["drug"]["drugindication"]
        #effect = row["patient"]["reaction"]["reactionmeddrapt"]
        term = row["term"]
        count = row["count"]
        #print term, count #indication, effect
        values.append((count, term))
    values = z_scorize_counts(values)
    values_eff = get_efficacy_values(values)
    return values, values_eff


def get_drugs_for_disease(disease):
    try:
        response = get_counts("disease-drug2", disease)
    except urllib2.HTTPError:
        print "No info for", disease
        return []
    values = []
    for row in response:
        term = row["term"]
        count = row["count"]
        values.append((count, term))
    #values = z_scorize_counts(values)
    return values


def get_diseases_for_drug(drug):
    values = []
    try:
        response = get_counts("drug-disease2", drug)
    except urllib2.HTTPError:
        print "No info for", drug
        return values
    for row in response:
        term = row["term"]
        count = row["count"]
        values.append((count, term))
    #values = z_scorize_counts(values)
    return values


if __name__ == "__main__":
    main()

