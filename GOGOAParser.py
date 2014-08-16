
import sys
import re

class GOGOAParser(object):
    """
    GO GOA Parser Class
    """
    name = "go_goa"
    description = "This program fills up tables in database biana related with gene ontology gene/protein annotations"

    db_name_to_biana_name = { "UniProtKB/Swiss-Prot": 'UniprotAccession', 
			    "UniProtKB/TrEMBL": 'UniprotAccession', 
			    "UniProtKB": 'UniprotAccession',
			    "ENSEMBL": 'ENSEMBL', 
			    "HINV": 'None', 
			    "TAIR": 'TAIR', 
			    "RefSeq": 'RefSeq', 
			    "VEGA": 'None',
			    "PDB": 'PDB',
			    'MGI': 'MGI'
			    }

    def __init__(self, filename):
	self.file_name = filename
	return

    def parse(self, exclude_evidences=None, id_type="genesymbol"):
	"""
	# EXP: Inferred from Experiment
	# IDA: Inferred from Direct Assay
	# IPI: Inferred from Physical Interaction
	# IMP: Inferred from Mutant Phenotype
	# IGI: Inferred from Genetic Interaction
	# IEP: Inferred from Expression Pattern 
	id_type: genesymbol | uniprot
	"""
	go_id_to_genes = {}
	if self.file_name.endswith(".gz") or self.file_name.endswith(".gz?rev=HEAD"):
	    import gzip
	    fd = gzip.open(self.file_name)
	else:
	    fd = open(self.file_name)
	for line in fd:
	    if line.startswith("!"):
		continue
	    words = line.rstrip('\n').split("\t")
	    if len(words) < 13:
		continue
	    db_name, db_id, gene_name, qualifier, go_id, evidence_ref, evidence_type, evidence_id, aspect, description, synonym, type, taxon_id = words[:13] #, date, curator, extension, form_id
    
	    if exclude_evidences is not None and evidence_type in exclude_evidences:
		continue

	    if any([ q=="NOT" for q in qualifier.split('|')]): # possible other qualifiers: colocalizes_with, contributes_to
		continue

	    if db_name not in self.db_name_to_biana_name:
		#print "Warning: xref db name is not recognized:", db_name
		db_name = None
	    else:
		db_name = self.db_name_to_biana_name[db_name]

	    if db_name is not None:
		xref = (db_name, db_id)

	    tax_ids = [ int(taxon.lstrip("taxon:")) for taxon in taxon_id.split('|') ]

	    if id_type == "uniprot":
		if db_name == "UniprotAccession":
		    for tax_id in tax_ids:
			go_id_to_genes.setdefault(go_id, set()).add((db_id, tax_id))
		else:
		    raise ValueError("Uknown db_name %s" % db_name) 
	    elif id_type == "genesymbol":
		for tax_id in tax_ids:
		    go_id_to_genes.setdefault(go_id, set()).add((gene_name, tax_id))
	    else:
		raise ValueError("Unknown id_type %s" % id_type)
	    
	fd.close()

	return go_id_to_genes


