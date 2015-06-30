import urllib2, os, cPickle, re, time
from bs4 import BeautifulSoup
from toolbox import text_utilities


def main():
    spl_id = "9c983da1-a8e0-479e-8816-5d5668ea242e" 
    #drug = "methotrexate"
    #disease = "type 2 diabetes mellitus"
    rxnorm_mapping_file = "/home/emre/arastirma/data/drug/fda/dailymed/rxnorm_mappings.txt"
    spl_id_to_names = get_rxnorm_mapping(rxnorm_mapping_file)
    print spl_id_to_names[spl_id]
    return 
    for spl in (spl_id,):
        response = get_data(parameter=spl)
        #print response
        name, indication, contraindication, warning = read_spl_data(response.read())
        print spl, name
        print indication
        print contraindication
        print warning
    return 


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
    for tag in soup.find_all('h1'):
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
    return name, indication, contraindication, warning

    
	
if __name__ == "__main__":
    main()

