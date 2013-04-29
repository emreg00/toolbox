#########################################################################
# PSI MI XML parser to parse phosphorylation info of proteins
#
# eg 15/04/2013
#########################################################################

from xml.etree.ElementTree import iterparse

def main():
    from time import clock
    parser = PSIMIXMLParser("/sbi/users/interchange/BIANA_FILES/biana_databases/biana_2011/NEW_BIANA_EXTERNAL_DATABASES/intact/human_small-27.xml") #../data/human_13.xml")
    t1 = clock()
    elements = parser.parse()
    t2 = clock()
    #print len(elements), elements[-1]
    print t2-t1 
    for element in elements:
	for interaction in element.interactions:
	    #print interaction.iid
	    participant_names=set()
	    id1, id2 = None, None
	    for participant in interaction.get_participants():
		#print element.interactors[participant].get_name()
		name = element.interactors[participant].get_name()
		if name == "TUBB":
		    id1 = participant
		elif name == "KRT81": #"PLOD2":
		    id2 = participant
		participant_names.add(name)
	    if "KRT81" in participant_names and "TUBB" in participant_names:
		print interaction.iid, id1, id2

    return

class PSIMIXMLParser(object):
    #NS="{net:sf:psidev:mi}" #"{http://psi.hupo.org/mi/mif}"

    def __init__(self, filename):
	self.file_name = filename
	#self.etree = ElementTree()
	return

    def parse(self):
	context = iterparse(self.file_name, ["start", "end"])
	context = iter(context)
	event, root = context.next()
	self.NS = "{%s}" % root.get("xmlns")
	stack = []
	stack.append(root)
	elements = []
	current_element = None
	current_interactor = None
	current_interaction = None
	for (event, elem) in context:
	    if event == "start":
		stack.append(elem.tag)
		if elem.tag == self.NS+"entry":
		    current_element = PSIMIXMLElement()
		elif elem.tag == self.NS+"interactorList":
		    if len(current_element.get_interactors()) > 0:
			print current_element.get_interactors()
		elif elem.tag == self.NS+"interactionList":
		    if len(current_element.get_interactions()) > 0:
			print current_element.get_interactions()
		elif elem.tag == self.NS+"interactor":
		    current_interactor = PSIMIXMLInteractor(elem.get("id"))
		elif elem.tag == self.NS+"interaction":
		    current_interaction = PSIMIXMLInteraction(elem.get("id"))
	    elif event == "end":
		top_tag = stack.pop()
		if elem.tag == self.NS+"interactor":
		    current_element.add_interactor(current_interactor)
		elif elem.tag == self.NS+"interaction":
		    current_element.add_interaction(current_interaction)
		elif elem.tag == self.NS+"alias" and stack[-2] == self.NS+"interactor":
		    if elem.get("type") == "gene name":
			current_interactor.add_name(elem.text)
		    elif elem.get("type") == "gene name synonym":
			current_interactor.add_synonym(elem.text)
		elif elem.tag == self.NS+"interactorRef":
		    current_interaction.add_participant(elem.text)
		elif elem.tag == self.NS+"entry":
		    elements.append(current_element) 
		    elem.clear()
	root.clear()
	return elements

class PSIMIXMLElement(object):

    def __init__(self):
	self.interactors = {}
	self.interactions = []

    def add_interactor(self, ielement):
	self.interactors[ielement.iid] = ielement

    def add_interaction(self, ielement):
	self.interactions.append(ielement)

    def get_interactors(self):
	return self.interactors

    def get_interactions(self):
	return self.interactions

    def __repr__(self):
	return "%s\t%s" % (self.interactors, self.interactions)

class PSIMIXMLInteractor(object):

    def __init__(self, iid):
	self.iid = iid
	self.name = None
	self.synonyms = []

    def add_name(self, name):
	self.name = name

    def add_synonym(self, synonym):
	self.synonyms.append(synonym)

    def get_name(self):
	return self.name

    def get_synonyms(self):
	return self.synonyms

    def __repr__(self):
	return "%s\t%s\t%s" % (self.iid, self.name, self.synonyms)

class PSIMIXMLInteraction(object):

    def __init__(self, iid):
	self.iid = iid
	self.participants = []

    def add_participant(self, participant):
	self.participants.append(participant)

    def get_participants(self):
	return self.participants

    def __repr__(self):
	return "%s\t%s" % (self.iid, self.participants)

if __name__ == "__main__":
    main()

