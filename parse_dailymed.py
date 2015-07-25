import urllib2, os, cPickle, re, time
from bs4 import BeautifulSoup
from toolbox import text_utilities, parse_drugbank


def main():
    spl_id = "e2158832-ad0b-4b9c-bbf0-1608a177bf85" #"29a6f213-bc93-4dbd-bab0-744722b6f0b8"
    spl_id2 = "9c983da1-a8e0-479e-8816-5d5668ea242e" 
    #drug = "methotrexate"
    #disease = "type 2 diabetes mellitus"
    rxnorm_mapping_file = "/home/emre/arastirma/data/drug/fda/dailymed/rxnorm_mappings.txt"
    spl_id_to_names = get_rxnorm_mapping(rxnorm_mapping_file)
    print spl_id_to_names[spl_id]
    #return 
    for spl in (spl_id, spl_id2):
	name, indication, contraindication, warning = fetch_spl_data("/home/emre/arastirma/data/drug/fda/dailymed/spl", spl)
        print spl, name
        print indication
        print contraindication
        print warning
    return 


def fetch_spl_data(output_dir, spl):
    out_file = output_dir + "/%s.html" % spl 
    response = get_data(parameter=spl)
    #print response
    if response is None:
	return None, None, None, None
    result = []
    f = open(out_file, 'w')
    for row in response:
	if row.strip() == "":
	    continue
    	f.write(row.replace("\r", ""))
    	result.append(row)
    f.close()
    name, indication, contraindication, warning = read_spl_data("".join(result)) # response.read()
    return name, indication, contraindication, warning


def get_drugbank_id_to_spl_id_mapping(rxnorm_mapping_file, name_to_drug, synonym_to_drug, selected_drugs=None, dump_file=None):
    if dump_file is not None and os.path.exists(dump_file):
	drugbank_id_to_spl_ids = cPickle.load(open(dump_file))
	return drugbank_id_to_spl_ids 
    spl_id_to_names = get_rxnorm_mapping(rxnorm_mapping_file)
    drugbank_id_to_spl_ids = {}
    for spl_id, names in spl_id_to_names.iteritems():
	drugbank_id = None
	for name in names:
	    if name.find(" / ") != -1: # Skip drug combinations
		continue
	    drugbank_id, drugbank_name = parse_drugbank.get_drugbank_id_from_name(name, name_to_drug, synonym_to_drug, regex_db_name = True) 
	    if drugbank_id is not None:
		break
	if drugbank_id is None:
	    continue
	#print name, drugbank_name
	drugbank_id_to_spl_ids.setdefault(drugbank_id, []).append(spl_id)
    cPickle.dump(drugbank_id_to_spl_ids, open(dump_file, 'w'))
    return drugbank_id_to_spl_ids


def get_rxnorm_mapping(mapping_file):
    """ 
    SETID|SPL_VERSION|RXCUI|RXSTRING|RXTTY
    """
    spl_id_to_names = {}
    f = open(mapping_file)
    header = f.readline()
    for line in f:
        words = line.strip().split("|")
        spl_id = words[0]
        name = words[3]
        spl_id_to_names.setdefault(spl_id, []).append(name)
    f.close()
    return spl_id_to_names


def get_data(command="setid", parameter=None):
    response = None
    url = 'http://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?%s=%s' % (command, parameter)
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


def read_spl_data(html_doc):
    name = None
    indication = [] 
    contraindication = [] 
    warning = []
    soup = BeautifulSoup(html_doc, "lxml")
    for tag in soup.find_all(['h1']):
        #print tag, tag.get("class")
        cls = tag.get("class")
        tag = tag.find_next()
        cls_inner = tag.get("class")  
        #print tag, cls
        if cls is not None and cls[0] == "Highlights":
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
        elif cls_inner is not None and cls_inner[0] == "long-title":
            if name is not None:
                print "Multiple name:", name
            name = tag.get_text().encode("ascii", "ignore")
            words = name.split(" - ")
            name = words[0].strip().lower()
    if len(indication) == 0:
	for tag in soup.find_all(['a']):
	    href = tag.get("href")
	    if href is not None and href[0] == "#":
		if tag.string is None:
		    continue
		header = tag.string.encode("ascii", "ignore").strip().lower()
		flag = None
		if header.find("contraindication") != -1:
		    flag = "contraindication"
		elif header.find("indication") != -1 or header in ("uses", "use", "usage"): # "indications", "indications and usage"
		    flag = "indication"
		if header.find("warning") != -1 and header.find("boxed") == -1:
		    flag = "warning"
		if flag is not None:
		    for tag_p in tag.find_all_next(["p", "li", "a"]):
			try:
			    if tag_p.name == "a":
				href = tag_p.get("href")
				if href is not None and href[0] == "#":
				    break
			except:
			    continue
			if tag_p.name == "p" or tag_p.name == "li":
			    txt = tag_p.get_text().strip()
			    if txt == "":
				continue
			    tag2 = tag_p.find_next("a")
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

