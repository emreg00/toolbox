###############################################################################
# To parse OBO formatted files
#
# Initial code contributed by Anand Patel from UCSD
###############################################################################

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
	elif tag_line.startswith("is_obsolete"):
	    self.type = "obsolete"
	    self.subset = tag_line.split(' ')[1]
        elif tag_line.startswith('alt_id'):
            parts = tag_line.split(' ')
	    self.type = parts[0].rstrip(':')
            self.id = parts[1]
        elif tag_line.startswith('xref'):
            parts = tag_line.split(' ')
	    self.type = parts[0].rstrip(':')
            self.id = parts[1]
	else:
	    self.type = "ignore"
        
    def __str__(self):
        return "relationship: %s %s ! %s"%(self.type, self.id, self.name)


def getOboGraph(fname, save_synonyms = False, exclude_obsolete = True):
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
            obo.node[id]['n'] = name
            # Alternative ids 
            obo.node[id]['alt_id'] = []
            # Cross references
            obo.node[id]['xref'] = []
            # Genes
            obo.node[id]['genes'] = []
            
            while True:
                stanza = obo_file.next().strip()
                if stanza=='':
                    break
		# Namespace 
		if stanza.startswith("namespace"):
		    namespace = stanza[11:]
		    obo.node[id]['t'] = namespace
		    continue
                stanza = Relationship(stanza, save_synonyms)
		if stanza.type == "ignore":
		    continue
		elif stanza.type == "subset":
		    if stanza.subset == "goslim_yeast":
			obo.node[id]['y'] = True
		    if stanza.subset.startswith("goslim"):
			obo.node[id]['a'] = True
		elif stanza.type == "obsolete":
		    if exclude_obsolete:
			continue
		    else:
			obo.node[id]['o'] = True
		# Synonyms 
		elif save_synonyms and stanza.type == "synonym":
		    obo.node[id].setdefault('s', set()).add(stanza.id)
		    #obo.node[id]['s'].add(stanza.id)
		elif stanza.type in ("alt_id", "xref"):
		    obo.node[id][stanza.type].append(stanza.id)
		elif stanza.type in ("is_a", "relationship"): 
		    obo.add_edge(id, stanza.id)
		    obo[id][stanza.id]['r'] = stanza.type
    
    return obo

