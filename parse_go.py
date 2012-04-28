from GO import GO

def main():
    species_prefix = "sgd" # "goa_human"
    go = GO("../data/gene_ontology.1_2-20_2_2012.obo",
	    False,
	    "../data/gene_association." + species_prefix + ".gz?rev=HEAD", #"../../../data/go/gene_association.goa_uniprot.gz",
	    exclude_evidences=[])
    print len(go.g.nodes()), len(go.go_id_to_genes)
    result = []
    for go_id, data in go.g.nodes(data=True):
	if go_id not in go.go_id_to_genes:
	    continue
	desc = data['n']
	result.append( (desc, [(gene, tax_id) for gene, tax_id in go.go_id_to_genes[go_id]] ) )
    return

if __name__ == "__main__":
    main()

