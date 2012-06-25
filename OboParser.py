
import networkx as nx
import re

class Relationship:
    def __init__(self, tag_line, save_synonyms = False):
        if tag_line.startswith('relationship'):
            #relationship: <type> <id> ! <name>
            parts = tag_line.split(' ')
            self.type = parts[1]
            self.id = parts[2]
            self.name = " ".join(parts[4:])
        elif tag_line.startswith('is_a'):
            parts = tag_line.split(' ')
	    self.type = parts[0].rstrip(':')
            self.id = parts[1]
            self.name = " ".join(parts[3:])
        elif save_synonyms and tag_line.startswith('synonym'):
	    temp = re.search("\"(.+)\"\s+(\w+)",tag_line)
	    #if temp.group(2) == "EXACT": #synonym = temp.group(1)
            self.type = "synonym" # synonym
            self.id = temp.group(1) # synonymous term
            self.name = temp.group(2) # synonym type
	elif tag_line.startswith("subset"):
	    self.type = "subset"
	    self.subset = tag_line.split(' ')[1]
	else:
	    self.type = "ignore"
        
    def __str__(self):
        return "relationship: %s %s ! %s"%(self.type, self.id, self.name)


def getOboGraph(fname, save_synonyms = False):
    obo_file = open(fname, 'rb')

    obo = nx.DiGraph()
    while True:
        try:
            line = obo_file.next().strip()
        except StopIteration:
            break
        if line=='[Term]':
            # Expect ID in the next line
            id = obo_file.next().strip()[4:]
            # Expect Name in the next line
            name = obo_file.next().strip()[6:]
            obo.add_node(id)
	    # Expect Namespace in the next line
            namespace = obo_file.next().strip()[11:]
            # Legible Name
            obo.node[id]['n'] = name
            obo.node[id]['t'] = namespace
            # Number of Arrays
            obo.node[id]['w'] = 0
            # Diff Express Genes
            obo.node[id]['g'] = {}
	    # Synonyms 
	    #obo.node[id]['s'] = set()
            
            while True:
                stanza = obo_file.next().strip()
                if stanza=='':
                    break
                
                stanza = Relationship(stanza, save_synonyms)
		#try:                
		if stanza.type == "ignore":
		    continue
		if stanza.type == "subset":
		    if stanza.subset == "goslim_yeast":
			obo.node[id]['y'] = True
		    continue
		if save_synonyms and stanza.type == "synonym":
		    obo.node[id].setdefault('s', set()).add(stanza.id)
		    #obo.node[id]['s'].add(stanza.id)
		    continue
		obo.add_edge(id, stanza.id)
		obo[id][stanza.id]['r'] = stanza.type
		#except:
		#    pass
    
    return obo


