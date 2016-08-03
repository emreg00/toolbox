#!/usr/bin/env python
# Cloned & modified from https://gist.github.com/yy/869845 
# see also original Python client from the webserver's authors
# at http://llama.mshri.on.ca/software.html
import sys
import json
import httplib

class FuncassociateClient(object):
    """Query funcassociate to find out enriched GO terms."""
    
    host = 'llama.mshri.on.ca'
    query_url = '/cgi/funcassociate/serv'

    def __init__(self):
        self.c = httplib.HTTPConnection(self.host)

    def close_conn(self):
        self.c.close()

    def jsonify(self, data):
        return (json.dumps(data)).encode('utf-8')

    def request(self, payload):
        self.c.request('POST', self.query_url, self.jsonify(payload), 
                       headers={'Content-type': 'application/json'})
        response = self.c.getresponse()
        if response.status == httplib.OK:
            return response.read()

    def available_species(self):
        payload = {'method': 'available_species', 
                   'id': 0,
                   }
        return self.request(payload)

    def available_namespaces(self, species=['Homo sapiens']):
        payload = {'method': 'available_namespaces',
                   'params': species,
                   'id': 123123123
                   }
        return self.request(payload)

    def go_associations(self, 
                        params=['Homo sapiens', 'uniprot_swissprot_accession', #'uniprot_id', 
                                ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
                               ]):
        payload = {'method': 'go_associations',
                   'params': params,
                   'id': 1,
                   }
        return self.request(payload)


    def functionate(self, query, species="Homo sapiens", namespace="entrezgene", genespace=None, mode="unordered", reps=2000, support=None, associations=None):
	params_dict = { "query": query, 
		       "species": species,
		       "namespace": namespace,
		       "mode": mode,
		       "reps": reps
			       }
	if associations is not None:
	    params_dict["support"] = support
	if associations is not None:
	    params_dict["associations"] = associations
	    del params_dict["species"]
	    del params_dict["namespace"]
	if genespace is not None:
	    params_dict["genespace"] = genespace
	payload = {'id': 1,
		   'method': 'functionate',
		   'params': [ params_dict ], 
		   'jsonrpc': '2.0'}
	response = self.request(payload=payload)
	#print response
        response = json.loads(response)
	#self.close_conn()
            
	if "error" in response and response["error"] is not None:
	    raise ValueError("Error %s" % response['error']['message'])

	response = response['result']
	return response


if __name__ == '__main__':
    functionate(["YLR406C", "YCR031C", "YDR471W"], "Saccharomyces cerevisiae", "sgd_systematic")
    fa = FuncassociateClient()
    print fa.available_species()
    print fa.available_namespaces()
    #print fa.go_associations()
    payload = {'id': 1,
               'method': 'functionate',
               'params': [ { "query": ["YLR406C", "YCR031C", "YDR471W"],
		   "species":"Saccharomyces cerevisiae",
		   "namespace": "sgd_systematic"
		   #"species":"Homo sapiens",
                           }], 
               'jsonrpc': '2.0'}
    print fa.request(payload=payload)
    fa.close_conn()
