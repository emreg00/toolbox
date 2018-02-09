import urllib2, json, urllib
import configuration
import ssl

CONFIG = configuration.Configuration() 

try:
    API_USER_KEY = CONFIG.get("OPENPHACTS_API_KEY") # change this value with custom API key
    API_APP_ID = CONFIG.get("OPENPHACTS_APP_ID")
except:
    print "Warning: OPENPHACTS API key not found!"
    API_USER_KEY = None
    API_APP_ID = None

FIELD_DRUG = "http://aers.data2semantics.org/resource/drug/"
FIELD_TARGET = "https://beta.openphacts.org/2.1/compound/pharmacology/pages?uri=" # proteinBinding
FIELD_COMPOUND = "https://beta.openphacts.org/2.1/mapUri?Uri=" # hasTarget


def main():
    drug = "paracetamol" 
    print get_drug_targets(drug)
    return


def get_data(command, parameter):
    if command == "target": 
        txt = '%s%s' % (urllib.quote_plus(FIELD_DRUG), parameter)
    else:
        raise ValueError("Unknown command: " + command)
    if API_USER_KEY is None:
        raise ValueError("API key missing!")
    else:
        url = FIELD_TARGET + '%s&app_id=%s&app_key=%s' % (txt, API_APP_ID, API_USER_KEY)
    #print url
    req = urllib2.Request(url)
    gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
    try:
	#response = urllib2.urlopen(req) # without context - also works
	response = urllib2.urlopen(req, context=gcontext)
    except urllib2.HTTPError:
	print "Problem with response (probably no info):", url 
        return None
    while True:
	try:
	    response = json.load(response)
	    break
	except:
	    print "Problem with response:", parameter
	    response = urllib2.urlopen(req) # without context
    return response["result"]


def get_drug_targets(drug):
    try:
        response = get_data("target", drug)
    except urllib2.HTTPError:
        print "No info for", drug
        return []
    values = []
    values_all = []
    for row in response["items"]:
	term = row["hasAssay"]["hasTarget"]
	#print term 
	if "targetOrganismName" in term and term["targetOrganismName"] == "Homo sapiens":
	    values.append(term["title"])
    return values, values_all


if __name__ == "__main__":
    main()

