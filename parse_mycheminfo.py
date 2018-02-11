import urllib2, json

base_url = "http://mychem.info/v1/"

def main():
    print get_resource_info()
    #drugbank_id = "DB00843"
    #drugbank_name = "donepezil"
    inchi_key = "ADEBPBSSDYVVLD-UHFFFAOYSA-N"
    print fetch_target_data(inchi_key, id_type="inchi_key", resource="drugbank", target_type="enzyme")
    print fetch_target_data(inchi_key, id_type="inchi_key", resource="drugcentral", target_type="enzyme")
    drugbank_id = "DB00316"
    drugbank_name = "Acetaminophen"
    #inchi_key = "RZVAJINKPMORJF-UHFFFAOYSA-N"
    print fetch_target_data(drugbank_name, id_type="db_name", resource="drugbank", target_type="target")
    print fetch_target_data(drugbank_id, id_type="db_id", resource="drugcentral", target_type="target")
    return 


def get_resource_info():
    response = get_data()
    response = json.load(response)
    return response["src_version"]


def fetch_target_data(drug, id_type="db_name", resource="drugbank", target_type="target"):
    """
    Return a list of triples containing uniprot id, action type and organism of 
	the targets of a given drug defined by drugbank name (db_name), 
	id (db_id) or inchi key (inchi_key)
    drug: descriptor of the drug in one of the allowed id types
    id_type: db_name | db_id | inchi_key
    resource: drugbank | drugcentral (currently the other resources returned in 
        get_resource_info are not supported)
    target_type: target (pharma action) | enzyme | transporter | carier | membrane receptor | ...
    """
    if id_type == "db_name":
	response = get_data("query", drug)
    elif id_type == "db_id" or id_type == "inchi_key":
	response = get_data("drug", drug)
    else:
	raise ValueError("Unrecognized id type: %s" % id_type)
    if response is None:
	return None
    response = json.load(response)
    targets_all = []
    if "hits" in response:
	for result in response["hits"]:
	    targets = parse_target_data_in_result(result, resource, target_type)
	    if targets is not None:
		targets_all.extend(targets)
    else:
	targets_all = parse_target_data_in_result(response, resource, target_type)
    return targets_all


def parse_target_data_in_result(response, resource, target_type):
    targets = []
    if resource not in response:
	return None
    response = response[resource]
    if resource == "drugbank":
	target_type = target_type + "s"
	if target_type not in response:
	    return targets
	for target in response[target_type]:
	    action = None
	    organism = None
	    if "uniprot" in target:
		try:
		    uniprot = target["uniprot"]
		except:
		    print response
		    print target
		if "actions" in target:
		    action = target["actions"]
		if "organism" in target:
		    organism = target["organism"].lower()
		targets.append((uniprot, action, organism))
    elif resource == "drugcentral":
	if "bioactivity" not in response:
	    return targets
	for target in response["bioactivity"]:
	    if target_type != "target" and target["target_class"].lower() != target_type:
		continue
	    action = None
	    organism = None
	    if "uniprot_id" in target:
		uniprot = target["uniprot_id"]
		if "action_type" in target:
		    action = target["action_type"].lower()
		if "swissprot" in target:
		    organism = target["swissprot"].split("_")[1].lower()
		targets.append((uniprot, action, organism))
    else:
	raise NotImplementedError()
    return targets


def get_data(command="metadata", parameter=None):
    """
    command: metadata | query | drug
    parameter: None (metadata) | drugbank name (query) | drugbank id or InChIKey (drug)
    """
    response = None
    if command == "metadata":
	url = base_url + '%s' % command
    elif command == "query":
	url = base_url + 'query?q=drugbank.name:%s' % parameter
    elif command == "drug":
	url = base_url + '%s/%s' % (command, parameter)
    print url
    #user_agent = 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'
    #values = {}
    #headers = { 'User-Agent' : user_agent }
    #data = urllib.urlencode(values)
    #req = urllib2.Request(url, "", headers)
    req = urllib2.Request(url)
    try:
        response = urllib2.urlopen(req)
    except:
        print "Problem with response:", parameter
    return response

	
if __name__ == "__main__":
    main()

