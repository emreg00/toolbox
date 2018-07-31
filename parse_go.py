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


def get_go_annotations(file_name, key_identifier="genesymbol", evidences_to_filter=None, tax_id=None):
    """
    key_identifier: genesymbol | uniprot
    tax_id: 9606
    evindeces_to_filter: ['IKR', 'HMP', 'HEP', 'EXP', 'TAS', 'ISS', 'IEP', 'IPI', 'ND', 'NAS', 'ISM', 'ISO', 'IGI', 'HDA', 'IBA', 'IEA', 'ISA', 'IC', 'IMP', 'IDA']
    """
    gene_to_go_ids = {}
    go_id_to_genes = {}
    evidences_to_filter = set()
    if evidences_to_filter is not None:
	evidences_to_filter = set(evidences_to_filter)
    molecules_all = set()
    evidences_all = set()
    with open(file_name) as f:
	for line in f.readlines():
	    if line[0] == "!":
		continue
	    words = line.strip("\n").split("\t")
	    if key_identifier == "genesymbol":
		gene = words[2]
	    else:
		gene = words[1]
	    go_term = words[4]
	    evidence = words[6]
	    genes = words[10]
	    molecule = words[11]
	    tax = words[12]
	    molecules_all.add(molecule)
	    evidences_all.add(evidence)
	    #print words
	    #print gene, go_term, evidence, molecule, tax, genes
	    if evidences_to_filter is not None and evidence in evidences_to_filter:
		continue
	    if tax != "taxon:"+tax_id:
		continue
	    go_id_to_genes.setdefault(go_term, set()).add(gene)
	    gene_to_go_ids.setdefault(gene, set()).add(go_term)
    #print molecules_all
    #print evidences_all
    return go_id_to_genes, gene_to_go_ids


if __name__ == "__main__":
    main()

