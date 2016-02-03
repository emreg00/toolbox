##############################################################################
# DrugBank v4.3 XML parser
#
# eg 15/01/2016
##############################################################################

from xml.etree.ElementTree import iterparse
import re, math

def main():
    base_dir = "/Users/eguney/Dropbox/sbnb/"
    drugbank_file = base_dir + "drugbank.xml" # test.xml
    get_drugs_for_targets(drugbank_file, "drugs.txt")
    return
    drugs_file = None #"drugs.txt"
    scores_file = "node_scores.sif.genesymbol"
    output_file = "drug_scores.txt"
    drug_to_targets, drug_to_description, drug_to_indication = get_drug_targets(drugbank_file, drugs_file)
    #print len(drug_to_targets)
    output_drug_targets(drug_to_targets)
    score_drugs_by_target_score(drug_to_targets, scores_file, output_file)
    return


class DrugBankXMLParser(object):
    NS="{http://drugbank.ca}"

    def __init__(self, filename):
	self.file_name = filename
	self.drug_to_name = {}
	self.drug_to_brands = {}
	self.drug_to_synonyms = {}
	self.drug_to_description = {}
	self.drug_to_groups = {}
	self.drug_to_indication = {}
	self.drug_to_interactions = {}
	self.drug_to_partner_ids = {}
	self.drug_to_partner_ids_active = {}
	self.drug_to_pubchem = {}
	self.drug_to_targets = {}
	self.drug_to_targets_active = {}
	self.partner_id_to_gene = {}
	self.partner_id_to_uniprot = {}
	return

    def parse(self, selected_names=None, exp=None):
	#exp1 = re.compile("(.+) [(].*[)]")
	#exp2 = re.compile("(.+) [\[].*[\]]")
	# get an iterable
	context = iterparse(self.file_name, ["start", "end"])
	# turn it into an iterator
	context = iter(context)
	# get the root element
	event, root = context.next()
	state_stack = [ root.tag ]
	drug = None
	drug_id = None
	drug_id_partner = None
	partner_id = None
	resource = None
	for (event, elem) in context:
	    if event == "start":
		state_stack.append(elem.tag)
		if elem.tag == self.NS+"partner":
		    if state_stack[-2] == self.NS+"partners":
			partner_id = elem.attrib["id"]
		elif elem.tag == self.NS+"resource":
		    uniprot_resource = False
		elif elem.tag == self.NS+"target":
		    if state_stack[-2] == self.NS+"targets":
			current_target = elem.attrib["partner"]
	    if event == "end":
		if elem.tag == self.NS+"name":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_name[drug_id] = elem.text
		if elem.tag == self.NS+"brand":
		    if state_stack[-2] == self.NS+"brands" and state_stack[-3] == self.NS+"drug":
			brand = elem.text
			#idx = brand.find("(")
			#if idx != -1:
			#    brand = brand[:idx]
			idx = brand.find(" [")
			if idx != -1:
			    brand = brand[:idx]
			brand = brand.strip().encode('ascii','ignore')
			if brand != "":
			    self.drug_to_brands.setdefault(drug_id, set()).add(brand) 
		if elem.tag == self.NS+"synonym":
		    if state_stack[-2] == self.NS+"synonyms" and state_stack[-3] == self.NS+"drug":
			synonym = elem.text
			#idx = synonym.find("(")
			#if idx != -1:
			#    synonym = synonym[:idx]
			idx = synonym.find(" [")
			if idx != -1:
			    synonym = synonym[:idx]
			synonym = synonym.strip().encode('ascii','ignore')
			if synonym != "":
			    self.drug_to_synonyms.setdefault(drug_id, set()).add(synonym) 
		if elem.tag == self.NS+"drugbank-id":
		    if state_stack[-2] == self.NS+"drug":
			if "primary" in elem.attrib:
			    drug_id = elem.text
		elif elem.tag == self.NS+"description":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_description[drug_id] = elem.text
		    if len(state_stack) > 3 and state_stack[-3] == self.NS+"drug-interactions" and state_stack[-2] == self.NS+"drug-interaction":
			self.drug_to_interactions[drug_id][drug_id_partner] = elem.text
		elif elem.tag == self.NS+"indication":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_indication[drug_id] = elem.text
		elif elem.tag == self.NS+"target":
		    if state_stack[-2] == self.NS+"targets":
			#current_target = elem.attrib["partner"]
			#! need to change
			self.drug_to_partner_ids.setdefault(drug_id, []).append(current_target)
		elif elem.tag == self.NS+"known-action":
		    if state_stack[-2] == self.NS+"target":
			if elem.text == "yes":
			    #! 
			    self.drug_to_partner_ids_active.setdefault(drug_id, set()).add(current_target)
		elif elem.tag == self.NS+"group":
		    if state_stack[-2] == self.NS+"groups":
			self.drug_to_groups.setdefault(drug_id, set()).add(elem.text)
		elif elem.tag == self.NS+"drug":
		    if len(state_stack) > 3 and state_stack[-3] == self.NS+"drug-interactions" and state_stack[-2] == self.NS+"drug-interaction":
			d = self.drug_to_interactions.setdefault(drug_id, {})
			drug_id_partner = elem.text
			d[drug_id_partner] = ""
		elif elem.tag == self.NS+"gene-name":
		    if state_stack[-3] == self.NS+"partners" and state_stack[-2] == self.NS+"partner":
			#!
			self.partner_id_to_gene[partner_id] = elem.text
		elif elem.tag == self.NS+"resource":
		    if state_stack[-3] == self.NS+"external-identifiers" and state_stack[-2] == self.NS+"external-identifier":
			resource = elem.text 
		elif elem.tag == self.NS+"identifier":
		    if state_stack[-3] == self.NS+"external-identifiers" and state_stack[-2] == self.NS+"external-identifier":
			if state_stack[-5] == self.NS+"partners" and state_stack[-4] == self.NS+"partner":
			    if resource == "UniProtKB":
				#! 
				self.partner_id_to_uniprot[partner_id] = elem.text
			elif state_stack[-4] == self.NS+"drug":
			    if resource == "PubChem Compound":
				self.drug_to_pubchem[drug_id] = elem.text
		elem.clear()
		state_stack.pop()
	root.clear()
        # Map target ids to uniprot ids
        for drug, partner_ids in self.drug_to_partner_ids.iteritems():
            for partner_id in partner_ids:
                try:
                    uniprot = self.partner_id_to_uniprot[partner_id]
                except:
                    # drug target has no uniprot
                    continue
                self.drug_to_targets.setdefault(drug, set()).add(uniprot)
		if drug in self.drug_to_partner_ids_active and partner_id in self.drug_to_partner_ids_active[drug]:
		    self.drug_to_targets_active.setdefault(drug, set()).add(uniprot)
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

