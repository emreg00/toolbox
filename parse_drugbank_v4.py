##############################################################################
# DrugBank v4.3 XML parser
#
# eg 15/01/2016
##############################################################################

from xml.etree.ElementTree import iterparse
import re, math

def main():
    #base_dir = "/Users/eguney/Dropbox/sbnb/"
    #file_name = base_dir + "drugbank.xml" 
    base_dir = "/home/eguney/data/drugbank/"
    file_name = base_dir + "test.xml"
    parser = DrugBankXMLParser(file_name)
    parser.parse()
    uniprot_to_drugs = {}
    for drug, target_to_action in parser.drug_to_target_to_action.iteritems():
	#print drug
	for uniprot in target_to_action:
	    #print uniprot
	    uniprot_to_drugs.setdefault(uniprot, set()).add(drug)
    print uniprot_to_drugs
    return


class DrugBankXMLParser(object):
    NS="{http://www.drugbank.ca}"

    def __init__(self, filename):
	self.file_name = filename
	self.drug_to_name = {}
	self.drug_to_description = {}
	self.drug_to_groups = {}
	self.drug_to_indication = {}
	self.drug_to_pharmacodynamics = {}
	self.drug_to_moa = {}
	self.drug_to_toxicity = {}
	self.drug_to_synonyms = {}
	self.drug_to_products = {}
	self.drug_to_brands = {}
	self.drug_to_uniprot = {}
	self.drug_to_interactions = {} 
	self.drug_to_pubchem = {}
	self.drug_to_pubchem_substance = {}
	self.drug_to_kegg = {}
	self.drug_to_kegg_compound = {}
	self.drug_to_pharmgkb = {}
	self.drug_to_target_to_action = {}
	self.drug_to_targets_paction = {} 
        self.drug_to_categories = {}
        self.drug_to_atc_codes = {}
        self.drug_to_inchi_key = {}
        self.drug_to_smiles = {}
	self.target_to_gene = {}
	self.target_to_uniprot = {}
	return

    def parse(self, selected_names=None, exp=None):
	# get an iterable
	context = iterparse(self.file_name, ["start", "end"])
	# turn it into an iterator
	context = iter(context)
	# get the root element
	event, root = context.next()
	state_stack = [ root.tag ]
	drug_id = None
	drug_id_partner = None
	current_target = None
	resource = None
        current_property = None 
	for (event, elem) in context:
	    if event == "start":
		state_stack.append(elem.tag)
		if elem.tag == self.NS+"drugbank-id":
		    if "primary" in elem.attrib and state_stack[-3] == self.NS+"drugbank" and state_stack[-2] == self.NS+"drug":
			drug_id = None
		    elif len(state_stack) > 3 and state_stack[-3] == self.NS+"drug-interactions" and state_stack[-2] == self.NS+"drug-interaction":
			drug_id_partner = None
		elif elem.tag == self.NS+"resource":
		    resource = None
		elif elem.tag == self.NS+"property":
                    current_property = None
		elif elem.tag == self.NS+"target": 
		    if state_stack[-2] == self.NS+"targets":
			current_target = None #elem.attrib["partner"]
	    if event == "end":
		if elem.tag == self.NS+"drugbank-id":
		    if state_stack[-2] == self.NS+"drug":
			if "primary" in elem.attrib:
			    drug_id = elem.text
			    #print drug_id
		    elif len(state_stack) > 3 and state_stack[-3] == self.NS+"drug-interactions" and state_stack[-2] == self.NS+"drug-interaction":
			d = self.drug_to_interactions.setdefault(drug_id, {})
			drug_id_partner = elem.text
			d[drug_id_partner] = ""
		elif elem.tag == self.NS+"name":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_name[drug_id] = elem.text.strip()
		    elif state_stack[-2] == self.NS+"product" and state_stack[-3] == self.NS+"products":
			product = elem.text
			product = product.strip().encode('ascii','ignore')
			if product != "":
			    self.drug_to_products.setdefault(drug_id, set()).add(product)
		    elif state_stack[-2] == self.NS+"international-brand" and state_stack[-3] == self.NS+"international-brands":
			brand = elem.text
			#idx = brand.find(" [")
			#if idx != -1:
			#    brand = brand[:idx]
			brand = brand.strip().encode('ascii','ignore')
			if brand != "":
			    self.drug_to_brands.setdefault(drug_id, set()).add(brand) 
		elif elem.tag == self.NS+"description":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_description[drug_id] = elem.text
		    if len(state_stack) > 3 and state_stack[-3] == self.NS+"drug-interactions" and state_stack[-2] == self.NS+"drug-interaction":
			self.drug_to_interactions[drug_id][drug_id_partner] = elem.text
		elif elem.tag == self.NS+"group":
		    if state_stack[-2] == self.NS+"groups":
			self.drug_to_groups.setdefault(drug_id, set()).add(elem.text)
		elif elem.tag == self.NS+"indication":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_indication[drug_id] = elem.text
		elif elem.tag == self.NS+"pharmacodynamics":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_pharmacodynamics[drug_id] = elem.text
		elif elem.tag == self.NS+"mechanism-of-action":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_moa[drug_id] = elem.text
		elif elem.tag == self.NS+"toxicity":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_toxicity[drug_id] = elem.text
		elif elem.tag == self.NS+"synonym":
		    if state_stack[-2] == self.NS+"synonyms" and state_stack[-3] == self.NS+"drug":
			synonym = elem.text
			idx = synonym.find(" [")
			if idx != -1:
			    synonym = synonym[:idx]
			synonym = synonym.strip().encode('ascii','ignore')
			if synonym != "":
			    self.drug_to_synonyms.setdefault(drug_id, set()).add(synonym) 
		elif elem.tag == self.NS+"category":
		    if state_stack[-2] == self.NS+"categories":
			self.drug_to_categories.setdefault(drug_id, set()).add(elem.text)
		elif elem.tag == self.NS+"atc-code":
		    if state_stack[-2] == self.NS+"atc-codes":
			self.drug_to_atc_codes.setdefault(drug_id, set()).add(elem.attrib["code"])
		elif elem.tag == self.NS+"id":	
		    if state_stack[-3] == self.NS+"targets" and state_stack[-2] == self.NS+"target":
			current_target = elem.text
			d = self.drug_to_target_to_action.setdefault(drug_id, {})
			d[current_target] = "unknown"
			#print current_target 
		elif elem.tag == self.NS+"action":	
		    if state_stack[-3] == self.NS+"target" and state_stack[-2] == self.NS+"actions":
			self.drug_to_target_to_action[drug_id][current_target] = elem.text
		elif elem.tag == self.NS+"known-action":
		    if state_stack[-2] == self.NS+"target":
			if elem.text == "yes":
			    self.drug_to_targets_paction.setdefault(drug_id, set()).add(current_target)
		elif elem.tag == self.NS+"gene-name":
		    if state_stack[-3] == self.NS+"targets" and state_stack[-2] == self.NS+"target":
			self.target_to_gene[current_target] = elem.text
		elif elem.tag == self.NS+"kind":
		    if state_stack[-3] == self.NS+"calculated-properties" and state_stack[-2] == self.NS+"property":
                        current_property = elem.text # InChIKey or SMILES
		elif elem.tag == self.NS+"value":
		    if state_stack[-3] == self.NS+"calculated-properties" and state_stack[-2] == self.NS+"property":
                        if current_property == "InChIKey":
                            inchi_key = elem.text # strip InChIKey=
                            if inchi_key.startswith("InChIKey="):
                                inchi_key = inchi_key[len("InChIKey="):]
                            self.drug_to_inchi_key[drug_id] = inchi_key
                        if current_property == "SMILES":
                            self.drug_to_smiles[drug_id] = elem.text 
		elif elem.tag == self.NS+"resource":
		    if state_stack[-3] == self.NS+"external-identifiers" and state_stack[-2] == self.NS+"external-identifier":
			resource = elem.text 
		elif elem.tag == self.NS+"identifier":
		    if state_stack[-3] == self.NS+"external-identifiers" and state_stack[-2] == self.NS+"external-identifier":
			if state_stack[-5] == self.NS+"target" and state_stack[-4] == self.NS+"polypeptide":
			    if resource == "UniProtKB":
				self.target_to_uniprot[current_target] = elem.text
			elif state_stack[-4] == self.NS+"drug":
			    if resource == "PubChem Compound":
				self.drug_to_pubchem[drug_id] = elem.text
			    elif resource == "PubChem Substance":
				self.drug_to_pubchem_substance[drug_id] = elem.text
			    elif resource == "KEGG Drug":
				self.drug_to_kegg[drug_id] = elem.text
			    elif resource == "KEGG Compound":
				self.drug_to_kegg_compound[drug_id] = elem.text
			    elif resource == "UniProtKB":
				self.drug_to_uniprot[drug_id] = elem.text
			    elif resource == "PharmGKB":
				self.drug_to_pharmgkb[drug_id] = elem.text
		elem.clear()
		state_stack.pop()
	root.clear()
	print self.drug_to_target_to_action #!
	print self.target_to_uniprot #!
        # Map target ids to uniprot ids
	drug_to_target_uniprots = {}
	drug_to_target_uniprots_paction = {}
        for drug, target_to_action in self.drug_to_target_to_action.iteritems():
            for target, action in target_to_action.iteritems():
                try:
                    uniprot = self.target_to_uniprot[target]
                except:
                    # drug target has no uniprot
                    continue
                drug_to_target_uniprots.setdefault(drug, set()).add(uniprot)
		if drug in self.drug_to_targets_paction and target in self.drug_to_targets_paction[drug]:
		    drug_to_target_uniprots_paction.setdefault(drug, set()).add(uniprot)
	self.drug_to_target_to_action = drug_to_target_uniprots
	self.drug_to_targets_paction = drug_to_target_uniprots_paction
	return 

    def get_synonyms(self, selected_drugs=None):
	name_to_drug = {}
	for drug, name in self.drug_to_name.iteritems():
	    if selected_drugs is not None and drug not in selected_drugs:
		continue
	    name_to_drug[name.lower()] = drug
	synonym_to_drug = {}
	for drug, synonyms in self.drug_to_synonyms.iteritems():
	    for synonym in synonyms:
		if selected_drugs is not None and drug not in selected_drugs:
		    continue
		synonym_to_drug[synonym.lower()] = drug
	for drug, brands in self.drug_to_brands.iteritems():
	    for brand in brands:
		if selected_drugs is not None and drug not in selected_drugs:
		    continue
		synonym_to_drug[brand.lower()] = drug
	return name_to_drug, synonym_to_drug


if __name__ == "__main__":
    main()

