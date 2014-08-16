import urllib2, os, cPickle, json

def main():
    #drug = "InChIKey=PSGAAPLEWMOORI-PEINSRQWSA-N" 
    drug = "InChIKey=WIGIZIANZCJQQY-UHFFFAOYSA-N" # Glimepiride
    cell_line = None #"MCF7"
    get_drug_signature(drug, cell_line)
    return


def get_data(command, parameter, parameter2=None):
    if command == "list":
        parameter = parameter.replace(" ", "+")
        txt = 'pertinfo?q={"inchi_key":"%s"}' % (parameter)
        # example: InChIKey=PSGAAPLEWMOORI-PEINSRQWSA-N
    elif command == "get":
        if parameter2 is None:
            txt = 'siginfo?q={"pert_id":"%s"}' % (parameter)
        else:
            txt = 'siginfo?q={"pert_id":"%s","cell_id":"%s"}' % (parameter, parameter2)
        # example: BRD-K82216340
    else:
        raise ValueError("Unknown command: " + command)
    url = 'http://api.lincscloud.org/a2/%s&user_key=lincsdemo' % txt 
    #print url
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    try:
        response = json.load(response)
    except:
        print "Problem with response:", parameter, parameter2
        response = None
    return response


def get_drug_signature(drug, cell_line):
    #try:
    response = get_data("list", drug)
    #except urllib2.HTTPError:
    if len(response) > 1:
        print "Warning: multiple perturbations"
        for row in response:
            pert_id = row["pert_id"]
            print pert_id
    elif len(response) == 0:
        print "No info for", drug
    	return None
    pert_id = response[0]["pert_id"]
    response = get_data("get", pert_id, cell_line)
    values = []
    perturbation_times = []
    for row in response:
        #print row["pert_time"], row["cell_id"]
        perturbation_times.append(int(row["pert_time"]))
    perturbation_time = max(perturbation_times)
    for row in response:
        if perturbation_time != int(row["pert_time"]):
            continue
        values_up = row["up100_full"] 
        values_down = row["dn100_full"] 
        #values_up = row["up50_lm"]
        #values_down = row["dn50_lm"]
        #print len(values_up), len(values_down)
        values.append((values_up, values_down))
    return values


def get_probe_mapping(probe_mapping_file, id_type="symbol"):
    """
	Parses trimmed HGU133A annotation file from GEO (GPL96)
	to get genes
	id_type: symbol / geneid
    """
    probe_to_genes = {}
    f = open(probe_mapping_file )
    line = f.readline()
    for line in f:
	probe, gene, geneid = line.strip("\n").split("\t")
	if id_type == "symbol":
	    for word in gene.split("///"):
		word = word.strip()
		probe_to_genes.setdefault(probe, []).append(word)
	elif id_type == "geneid":
	     for word in geneid.split("///"):
		word = word.strip()
		probe_to_genes.setdefault(probe, []).append(word)
	else:
	    raise ValueError("Unknown id_type " + id_type)
    f.close()
    return probe_to_genes


if __name__ == "__main__":
    main()

