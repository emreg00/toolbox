import urllib2, os, cPickle, json

def main():
    drug = "InChIKey=PSGAAPLEWMOORI-PEINSRQWSA-N" 
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
    response = json.load(response)
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
    for row in response:
        print row["pert_time"], row["cell_id"]
        values_up = row["up100_full"] #"up50_lm"]
        values_down = row["dn100_full"] #"dn50_lm"]
        print len(values_up), len(values_down)
    return 


if __name__ == "__main__":
    main()

