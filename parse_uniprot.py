#########################################################################
# Uniprot XML parser to parse phosphorylation info of proteins
#
# eg 29/07/2009
#########################################################################

#from xml.etree.ElementTree import ElementTree
from xml.etree.ElementTree import iterparse
import TsvReader

def main():
    file_name = "../data/disease/uniprot/humdisease.txt"
    mim_to_mesh_values = get_mim_to_mesh(file_name)
    print len(mim_to_mesh)
    print mim_to_mesh["600807"]
    return
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


def get_uniprot_to_geneid(file_name, uniprot_ids=None, only_min=True):
    """
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
    zcat HUMAN_9606_idmapping_selected.tab.gz | cut -f 1,3 > idmapping.tab
    To parse HUMAN_9606_idmapping.dat file (trimmed to two columns) from Uniprot 
    """
    uniprot_to_geneid = {}
    #geneid_to_uniprots = {}
    f = open(file_name)
    f.readline()
    for line in f:
	uniprot, geneid = line.split("\t")
	geneid = geneid.strip()
	uniprot = uniprot.strip()
	if geneid == "" or uniprot == "":
	    continue
	if uniprot_ids is not None and uniprot not in uniprot_ids:
	    continue
	if only_min:
	    geneid = str(min(map(int, geneid.split("; "))))
	uniprot_to_geneid[uniprot] = geneid
    f.close()
    return uniprot_to_geneid


def get_uniprot_to_geneid_from_idmapping_file(file_name, uniprot_ids=None):
    """
    To parse idmapping.tab from Uniprot 
    Useful for id mapping of non-human species
    """
    parser = TsvReader.TsvReader(file_name, delim="\t", inner_delim=";")
    column_to_index, id_to_values = parser.read(fields_to_include=["UniProtKB-AC", "GeneID (EntrezGene)"], keys_to_include=uniprot_ids, merge_inner_values=True)
    uniprot_to_geneid = {}
    for uniprot, values in id_to_values.iteritems():
    	for val in values:
    	    geneid = val[column_to_index["geneid (entrezgene)"]]
	    #if uniprot in uniprot_to_geneid:
	    #	print "multiple gene id", uniprot
	    #uniprot_to_geneid.setdefault(uniprot, set()).add(geneid)
	    uniprot_to_geneid[uniprot] = geneid
    return uniprot_to_geneid


def get_mim_to_mesh(file_name):
    """
    To parse humdisease.txt from Uniprot
    """
    mim_to_mesh_values = {}
    f = open(file_name)
    line = f.readline()
    while not line.startswith("ID"):
	line = f.readline()
    words = line.strip().split()
    disease = " ".join(words[1:]).rstrip(".")
    for line in f:
	words = line.strip().split()
	if words[0] == "ID":
	    disease = " ".join(words[1:]).rstrip(".")
	if words[0] == "DR":
	    id_type = words[1].lower().rstrip(";")
	    if id_type == "mesh":
		mesh = words[2].rstrip(".")
	    elif id_type == "mim":
		mim = words[2].rstrip(";")
	if line.startswith("//"):
	    #if mim in mim_to_mesh_values and mim_to_mesh_values[mim][1] == mesh:
		#continue
	    #if mim in mim_to_mesh_values: print mim, mim_to_mesh_values[mim], disease, mesh
	    mim_to_mesh_values.setdefault(mim, []).append((disease, mesh))
    f.close()
    return mim_to_mesh_values

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

