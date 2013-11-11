##############################################################################
# DrugBank XML parser to parse targets of drugs
#
# eg 10/11/2011
##############################################################################

from xml.etree.ElementTree import iterparse
import re, math

def main():
    drugbank_file = "../data/drugbank.xml" # test.xml
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
	self.drug_to_pubchem = {}
	self.drug_to_targets = {}
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
		if elem.tag == self.NS+"resource":
		    uniprot_resource = False
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
			self.drug_to_partner_ids.setdefault(drug_id, []).append(elem.attrib["partner"])
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
			self.partner_id_to_gene[partner_id] = elem.text
		elif elem.tag == self.NS+"resource":
		    if state_stack[-3] == self.NS+"external-identifiers" and state_stack[-2] == self.NS+"external-identifier":
			resource = elem.text 
		elif elem.tag == self.NS+"identifier":
		    if state_stack[-3] == self.NS+"external-identifiers" and state_stack[-2] == self.NS+"external-identifier":
			if state_stack[-5] == self.NS+"partners" and state_stack[-4] == self.NS+"partner":
			    if resource == "UniProtKB":
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
	return 

def get_drugs_for_targets(file_name, output_file):
    parser = DrugBankXMLParser(file_name)
    parser.parse()
    uniprot_to_drugs = {}
    for drug, targets in parser.drug_to_targets.iteritems():
	#print drug
	for uniprot in targets:
	    uniprot_to_drugs.setdefault(uniprot, set()).add(drug)
    f = open(output_file, 'w')
    for uniprot, drugs in uniprot_to_drugs.iteritems():
	f.write("%s\t%s\n" % (uniprot, ";".join(drugs)))
    f.close()
    return

def output_drug_info(file_name, output_file):
    parser = DrugBankXMLParser(file_name)
    parser.parse()
    f = open(output_file, 'w')
    f.write("drugbank id\tname\tgroups\tpubchem id\tdescription\tindication\ttargets\n")
    for drug, name in parser.drug_to_name.iteritems():
        name = name.encode('ascii','ignore')
        try:
            groups = parser.drug_to_groups[drug]
        except:
            groups = []
        try:
            description = parser.drug_to_description[drug]
            description = description.replace("\n", "").encode('ascii','ignore')
        except:
            description = ""
        try:
            indication = parser.drug_to_indication[drug]
            indication = indication.replace("\n", "").encode('ascii','ignore')
        except:
            #print drug
            indication = ""
	if drug in parser.drug_to_pubchem:
            pubchem = parser.drug_to_pubchem[drug]
        else:
	    pubchem = "" 
        if drug in parser.drug_to_targets:
            targets = parser.drug_to_targets[drug]
        else:
            targets = []
        try:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, name, ";".join(groups), pubchem, description, indication, ";".join(targets)))
        except:
            print drug, name, groups, pubchem, description, indication, targets
            return
    f.close()
    return

def get_drug_info(drug_info_file):
    drug_to_values = {}
    f = open(drug_info_file)
    header = f.readline().strip().split("\t")
    col_to_idx = dict((k, i) for i, k in enumerate(header[1:]))
    for line in f:
        words = line.strip("\n").split("\t")
        drug_to_values[words[0]] = words[1:]
    return col_to_idx, drug_to_values

def get_disease_specific_drugs(parser, phenotypes):
    import re
    #keywords = ["cancer", "neoplasm", "leukemia", "lymphoma", "carcinoma", "melanoma"]
    #exps = [ re.compile(keyword) for keyword in keywords ]
    exps = [ re.compile(keyword.lower()) for keyword in phenotypes ]
    disease_to_drugs = {}
    for drug, indication in parser.drug_to_indication.iteritems():
	if indication is None:
	    continue
	#if any(map(lambda x: x is not None, [ exp.search(indication) for exp in exps ])):
	    #disease = keywords[0]
	    #disease_to_drugs.setdefault(disease, set()).add(drug)
	for disease, exp in zip(phenotypes, exps):
	    if exp.search(indication) is not None:
		disease_to_drugs.setdefault(disease, set()).add(drug)
    return disease_to_drugs


def get_drug_targets(file_name, drugs_file=None):
    parser = DrugBankXMLParser(file_name)
    parser.parse()
    drugs = None
    if drugs_file is not None:
	drugs = set([ line.strip().lower() for line in open(drugs_file) ])
	exp = re.compile("brain")
	exp2 = re.compile("metastasis")
	for drug, description in parser.drug_to_description.iteritems():
	    #drug = drug.lower()
	    if description is None:
		continue
	    m = exp.search(description)
	    m2 = exp2.search(description)
	    if True: #! m is not None and m2 is not None:
		drugs.add(drug)
	for drug, indication in parser.drug_to_indication.iteritems():
	    #drug = drug.lower()
	    if indication is None:
		continue
	    m = exp.search(indication)
	    m2 = exp2.search(indication)
	    if True: #! m is not None and m2 is not None:
		drugs.add(drug)
	#print drugs

    drug_to_targets = {}
    for drug, partner_ids in parser.drug_to_partner_ids.iteritems():
	#drug = drug.lower()
	if drugs is not None and drug not in drugs: 
	    continue
	#print drug
	for partner_id in partner_ids:
	    gene = parser.partner_id_to_gene[partner_id]
	    if gene is None:
		continue
	    drug_to_targets.setdefault(drug, set()).add(gene)
    return drug_to_targets, parser.drug_to_description, parser.drug_to_indication

def output_drug_targets(drug_to_targets):
    f = open("drug_to_targets.txt", 'w')
    f2 = open("drug_targets.txt", 'w')
    for drug, targets in drug_to_targets.iteritems():
	f.write("%s\t%s\n" % (drug, "\t".join(targets)))
	f2.write("%s\n" % "\n".join(targets))
    f.close()
    f2.close()
    return

def score_drugs_by_target_score(drug_to_targets, scores_file, output_file):
    gene_to_score = dict([ line.strip().split() for line in open(scores_file)])
    values = []
    for drug, targets in drug_to_targets.iteritems():
	scores = []
	for target in targets:
	    if target in gene_to_score:
		scores.append(float(gene_to_score[target]))
	if len(scores) == 0:
	    continue
	values.append((calculate_score(scores), drug))
    values.sort()
    values.reverse()
    f = open(output_file, 'w')
    for score, drug in values:
	f.write("%s\t%s\n" % (drug, str(score)))
    f.close()
    return 

def calculate_drug_score_from_targets(values):
    val = 0.0
    for value in values:
	val += value * value
    return math.sqrt(val)

if __name__ == "__main__":
    main()

