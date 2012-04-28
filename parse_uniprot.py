#########################################################################
# Uniprot XML parser to parse phosphorylation info of proteins
#
# eg 29/07/2009
#########################################################################

#from xml.etree.ElementTree import ElementTree
from xml.etree.ElementTree import iterparse

def main():
    from time import clock
    parser = UniprotXMLParser("../data/Q12888.xml")
    #parser = UniprotXMLParser("../../data/phosphorylation/uniprot/uniprot-phosphorylation-large-scale-analysis.xml")
    #ids = parser.parse_ids()
    #print map(len, ids)
    #print ids[-1]
    t1 = clock()
    elements = parser.parse()
    t2 = clock()
    print len(elements), elements[-1]
    print t2-t1 
    return

class UniprotXMLParser(object):
    NS="{http://uniprot.org/uniprot}"
    psiteDesc_to_psiteChar = {  "Phosphoserine": "S",
				"Phosphothreonine": "T",
				"Phosphotyrosine": "Y",
				"Phosphohistidine": "H" }

    def __init__(self, filename):
	self.file_name = filename
	#self.etree = ElementTree()
	return

    def parse_ids_high_mem(self):
	self.etree = ElementTree()
	tree = self.etree.parse(self.file_name)
	#ids = tree.findall(self.NS+"accession")
	ids = []
	sub_ids = None
	for e in tree.getiterator():
	    if e.tag == self.NS+"entry":
		if sub_ids is not None:
		    ids.append(sub_ids)
		sub_ids = []
	    if e.tag == self.NS+"accession":
		sub_ids.append(e.text)
	ids.append(sub_ids)
	return ids

    def parse_ids(self):
	ids = []
	sub_ids = []
	# get an iterable
	context = iterparse(self.file_name, ["start", "end"])
	# turn it into an iterator
	context = iter(context)
	# get the root element
	event, root = context.next()
	for (event, elem) in context:
	    if event == "end":
		if elem.tag == self.NS+"accession":
		    sub_ids.append(elem.text)
		if elem.tag == self.NS+"entry":
		    ids.append(sub_ids)
		    sub_ids = []
		    elem.clear()
	root.clear()
	return ids

    def parse(self):
	ignored_modification_types = set()
	context = iterparse(self.file_name, ["start", "end"])
	context = iter(context)
	event, root = context.next()
	elements = []
	current_element = None
	current_position = None
	for (event, elem) in context:
	    if event == "start":
		if elem.tag == self.NS+"entry":
		    current_element = UniprotXMLElement()
	    elif event == "end":
		if elem.tag == self.NS+"accession":
		    current_element.add_id(elem.text)
		elif elem.tag == self.NS+"organism":
		    db_elm = elem.find(self.NS+"dbReference") #only looks at sublevel - alternative: keep tag stack
		    if db_elm.get("type") == "NCBI Taxonomy":
			current_element.set_tax(db_elm.get("id"))
		elif elem.tag == self.NS+"feature" and elem.get("type") == "modified residue":
		    #print elem.getchildren()
		    #pos_elm = elem.find(self.NS+"position")
		    #if elem.get("status") == "probable":
		    #	continue
		    for sub_elm in elem.getiterator():
			if sub_elm.tag == self.NS+"position":
			    pos_elm = sub_elm
		    pos = pos_elm.get("position")
		    desc = elem.get("description")
		    vals = desc.split(";")
		    type = vals[0]
		    kinase = vals[1][vals[1].find("by")+2:].strip() if (len(vals) > 1) else None
		    if self.psiteDesc_to_psiteChar.has_key(type):
			type = self.psiteDesc_to_psiteChar[type]
			current_element.add_psite(pos, type, kinase)
		    else:
			ignored_modification_types.add(type)
		elif elem.tag == self.NS+"entry":
		    seq_elm = elem.find(self.NS+"sequence")
		    current_element.set_sequence(seq_elm.text)
		    elements.append(current_element) 
		    elem.clear()
	root.clear()
	print "Ignored mofications: ", ignored_modification_types
	return elements

class UniprotXMLElement(object):

    def __init__(self):
	self.ids = []
	self.taxid = None
	self.phosphosites = []
	self.sequence = None

    def add_id(self, id):
	self.ids.append(id)

    def set_tax(self, taxid):
	self.taxid = taxid

    def add_psite(self, pos, type=None, kinase=None):
	self.phosphosites.append( (pos, type, kinase) )

    def set_sequence(self, seq):
	self.sequence = seq.replace("\n","")

    def get_ids(self):
	return self.ids

    def get_tax(self):
	return self.taxid

    def get_psites(self):
	return self.phosphosites

    def get_sequence(self):
	return self.sequence

    def __repr__(self):
	return "%s\t%s\t%s\t%s" % (self.ids, self.taxid, self.phosphosites, self.sequence)


if __name__ == "__main__":
    main()

