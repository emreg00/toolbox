##############################################################################
# DrugBank XML parser to parse targets of drugs
#
# eg 10/11/2011
##############################################################################

from xml.etree.ElementTree import iterparse
import re, math

def main():
    drugbank_file = "../data/drugbank.xml" # test.xml
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
	self.drug_to_description = {}
	self.drug_to_indication = {}
	self.drug_to_partner_ids = {}
	self.partner_id_to_gene = {}
	return

    def parse(self, selected_names=None, exp=None):
	# get an iterable
	context = iterparse(self.file_name, ["start", "end"])
	# turn it into an iterator
	context = iter(context)
	# get the root element
	event, root = context.next()
	state_stack = [ root.tag ]
	drug = None
	partner_id = None
	for (event, elem) in context:
	    if event == "start":
		state_stack.append(elem.tag)
		if elem.tag == self.NS+"partner":
		    if state_stack[-2] == self.NS+"partners":
			partner_id = elem.attrib["id"]
	    if event == "end":
		if elem.tag == self.NS+"name":
		    if state_stack[-2] == self.NS+"drug":
			drug = elem.text
		elif elem.tag == self.NS+"description":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_description[drug] = elem.text
		elif elem.tag == self.NS+"indication":
		    if state_stack[-2] == self.NS+"drug":
			self.drug_to_indication[drug] = elem.text
		elif elem.tag == self.NS+"target":
		    if state_stack[-2] == self.NS+"targets":
			self.drug_to_partner_ids.setdefault(drug, []).append(elem.attrib["partner"])
		elif elem.tag == self.NS+"gene-name":
		    if state_stack[-3] == self.NS+"partners" and state_stack[-2] == self.NS+"partner":
			self.partner_id_to_gene[partner_id] = elem.text
		elem.clear()
		state_stack.pop()
	root.clear()
	return 

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

def calculate_score(values):
    val = 0.0
    for value in values:
	val += value * value
    return math.sqrt(val)

if __name__ == "__main__":
    main()

