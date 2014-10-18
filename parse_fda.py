import urllib2, os, cPickle, json
from stat_utilities import calc_mean_and_sigma
import configuration

CONFIG = configuration.Configuration() 

API_USER_KEY = CONFIG.get("FDA_API_KEY") #None # change this value with custom API key
LIMIT = 100
FIELD_DRUG = "patient.drug.medicinalproduct"
FIELD_DISEASE = "patient.drug.drugindication"
FIELD_EFFECT = "patient.reaction.reactionmeddrapt"

def main():
    #drug = "montelukast" 
    #drug = "vitamin e" # "vit. e"
    #drug = "donepezil"
    drug = "methotrexate"
    #disease = "type 2 diabetes"
    #disease = "alzheimer"
    #disease = "asthma"
    #disease = "acute lymphocytic leukaemia"
    disease = "rheumatoid arthritis"
    #disease = "bone sarcoma"
    condition = None #"drug ineffective"
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


def get_data_helper(command, parameter, parameter2=None, parameter_effect=None, skip=0):
    parameter = parameter.replace(" ", "+")
    if command == "drug":
        txt = '%s:"%s"' % (FIELD_DRUG, parameter)
    elif command == "disease":
        txt = '%s:"%s"' % (FIELD_DISEASE, parameter)
    elif command == "drug-disease":
	assert parameter2 is not None
	parameter2 = parameter2.replace(" ", "+")
	txt = '%s:"%s"+AND+%s:"%s"' % (FIELD_DRUG, parameter, FIELD_DISEASE, parameter2)
    elif command == "drug-disease-effect":
	assert (parameter2 is not None and parameter_effect is not None)
	parameter2 = parameter2.replace(" ", "+")
	parameter_effect = parameter_effect.replace(" ", "+")
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
    parameter = parameter.replace(" ", "+")
    if parameter2 is not None:
        parameter2 = parameter2.replace(" ", "+")
    if parameter_effect is not None:
        parameter_effect = parameter_effect.replace(" ", "+")
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
	print "No info for", parameter, parameter2, parameter_effect
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
    val = parameter.lower().replace("+", " ")
    if command in ("disease-drug", "disease-drug-effect"):
	val = parameter2.lower().replace("+", " ")
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
    m, s = calc_mean_and_sigma(zip(*count_term_pairs)[0])
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
    try:
        response = get_counts("drug-disease2", drug)
    except urllib2.HTTPError:
        print "No info for", drug
        return []
    values = []
    for row in response:
        term = row["term"]
        count = row["count"]
        values.append((count, term))
    #values = z_scorize_counts(values)
    return values


if __name__ == "__main__":
    main()

