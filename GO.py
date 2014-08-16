
from OBO import OBO

class GO(OBO):

    def __init__(self, file_name, save_synonyms = False, go_goa_file = None, exclude_evidences=None, id_type="genesymbol"):
	OBO.__init__(self, file_name, save_synonyms)
	self.go_id_to_genes = None
        self.gene_to_go_ids = None
	if go_goa_file is not None:
	    self._get_classification(go_goa_file, exclude_evidences, id_type)
	self.tf_genes = None
	self.tf_related_genes = None
	return

    def get_classification(self):
        """
            go_id_to_genes: go_id => (gene, tax_id)
        """
	if self.go_id_to_genes is None:
	    raise ValueError("GOA file not provided during initialization")
	return self.go_id_to_genes
	
    def _get_classification(self, go_goa_file, exclude_evidences, id_type):
        """
            go_id_to_genes: go_id => (gene, tax_id)
        """
	from GOGOAParser import GOGOAParser
	parser = GOGOAParser(go_goa_file)
	self.go_id_to_genes = parser.parse(exclude_evidences, id_type)
	return 

    def get_tf_genes(self):
	"""
	GO:0003700 transcription factor activity
	"""
	self.get_ontology_extended_id_mapping()
	if self.tf_genes is None:
	    self.tf_genes = self.get_genes("GO:0003700")
	return self.tf_genes

    def get_tf_related_genes(self):
	"""
	GO:0003700 : transcription factor activity [969 gene products]
	GO:0030528 : transcription regulator activity [1528 gene products]
	GO:0045449 : regulation of transcription [2654 gene products]
	GO:0008134 : transcription factor binding [524 gene products]
	"""
	self.get_ontology_extended_id_mapping()
	if self.tf_related_genes is None:
	    self.tf_related_genes = self.get_genes("GO:0003700")
	    self.tf_related_genes |= self.get_genes("GO:0030528")
	    self.tf_related_genes |= self.get_genes("GO:0045449")
	    self.tf_related_genes |= self.get_genes("GO:0008134")
	return self.tf_related_genes


    def get_genes(self, go_id, include_descendants = True):
	go_ids = set([go_id])
	if include_descendants:
	    self.get_ontology_extended_id_mapping()
	    if self.child_to_parent.has_key(go_id):
		[ go_ids.add(go_sub) for go_sub in
		    self.child_to_parent[go_id] ]
	go_genes = set()
	for go_sub in go_ids:
	    if self.go_id_to_genes.has_key(go_sub):
		go_genes |= self.go_id_to_genes[go_sub]
	return go_genes

    def get_go_terms_of_gene(self, gene, namespace=None):
        """
        namespace: None | biological_process | molecular_function | cellular_compartment
        """
        if self.gene_to_go_ids is None:
            go_id_to_genes = self.get_classification()
            self.gene_to_go_ids = {}
            for go_id, gene_and_taxids in go_id_to_genes.iteritems():
                if namespace is not None:
                    if go_id not in self.g.node or self.g.node[go_id]['t'] != namespace: # go_id may be in alt_id
                        continue
                for gene_symbol, tax_id in gene_and_taxids:
                    self.gene_to_go_ids.setdefault(gene_symbol, set()).add(go_id)
        go_ids = set()
        if gene in self.gene_to_go_ids:
            go_ids = self.gene_to_go_ids[gene]
        return go_ids


def main(go_fname='gene_ontology.1_2.obo',
         go_goa_fname='gene_association.goa_human',
         exclude_evidences=[]):
    go = GO(go_fname, False, go_goa_fname, 
            exclude_evidences=exclude_evidences)
    #print go.get_tf_genes()
    return go

if __name__=='__main__':
    main()
