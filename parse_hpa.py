
import parse_ncbi

def main():
    mapping_file = "/home/emre/data/ncbi/geneid_to_symbol.txt"
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    file_name = "/home/emre/data/tissue/hpa_subcellular_location.tsv"
    geneid_to_localization = get_geneid_to_localization(file_name, mapping_file) # , only_extracellular=True
    print geneid_to_localization["8648"] # NCOA1
    print geneid_to_localization["207"] # AKT1
    file_name_membrane = "/home/emre/data/tissue/hpa_predicted_membrane.tsv"
    file_name_secreted = "/home/emre/data/tissue/hpa_predicted_secreted.tsv"
    geneid_to_localization = get_geneid_to_predicted_localization(file_name_membrane, file_name_secreted, mapping_file)
    print geneid_to_localization["5925"] # RB1
    return


def get_geneid_to_localization(file_name, mapping_file, only_extracellular=False, filter_uncertain=False):
    """
    Relevant extracellular term: Plasma membrane (GO:0005886)
    Some other terms: Vesicles (GO:0043231);Nuclear bodies (GO:0016604);Nucleoplasm (GO:0005654);Cytosol (GO:0005829);Nucleoli (GO:0005730);Mitochondria (GO:0005739);Cell Junctions (GO:0030054);Endoplasmic reticulum (GO:0005783);Centrosome (GO:0005813);Nucleoli fibrillar center (GO:0001650);Golgi apparatus (GO:0005794);Nuclear speckles (GO:0016607)
    Reliability types: Uncertain | Approved | Supported | Enhanced
    """
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    geneid_to_localization = {}
    ids_unmatched = set()
    f = open(file_name)
    f.readline()
    for line in f:
	words = line.strip("\n").split("\t")
	# Evidence and reliability are for experimental characterization
	gene, reliability, go_terms = words[1], words[2], words[-1]
	if filter_uncertain and reliability == "Uncertain":
	    continue
	if gene in name_to_geneid:
	    geneid = name_to_geneid[gene]
	else:
	    #print "Unmatched id", gene
	    ids_unmatched.add(gene)
	    continue
	for go in go_terms.split(";"):
	    idx = go.find("(")
	    go = go[idx+1:-1]
	    if only_extracellular and go != "GO:0005886":
		continue
	    geneid_to_localization.setdefault(geneid, set()).add((go, reliability))
    f.close()
    print "Unmatched:", len(ids_unmatched) #, ", ".join(sorted(ids_unmatched))
    return geneid_to_localization


def get_geneid_to_predicted_localization(file_name_membrane, file_name_secreted, mapping_file):
    """
    Relevant protein classes: Predicted secreted proteins | Predicted membrane proteins
    Note that evidence and reliability are for experimental characterization
    """
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    geneid_to_localization = {}
    ids_unmatched = set()
    f = open(file_name_membrane)
    f.readline()
    for line in f:
	words = line.strip("\n").split("\t")
	gene, protein_class, evidence, reliability = words[0], words[6], words[7], words[11]
	#if reliability not in ("Approved", "Supported", "Enhanced"):
	#    continue
	if gene in name_to_geneid:
	    geneid = name_to_geneid[gene]
	else:
	    #print "Unmatched id", gene
	    ids_unmatched.add(gene)
	    continue
	for localization in protein_class.split(", "):
	    localization = localization.replace("Predicted ", "").replace(" proteins", "").replace(" genes", "")
	    if localization == "membrane":
		geneid_to_localization.setdefault(geneid, set()).add((localization, reliability))
    f.close()
    f = open(file_name_secreted)
    f.readline()
    for line in f:
	words = line.strip("\n").split("\t")
	# Evidence and reliability are for experimental characterization
	gene, protein_class, evidence, reliability = words[0], words[6], words[7], words[11]
	if gene in name_to_geneid:
	    geneid = name_to_geneid[gene]
	else:
	    #print "Unmatched id", gene
	    ids_unmatched.add(gene)
	    continue
	for localization in protein_class.split(", "):
	    localization = localization.replace("Predicted ", "").replace(" proteins", "").replace(" genes", "")
	    if localization == "secreted":
		geneid_to_localization.setdefault(geneid, set()).add((localization, reliability))
    f.close()
    print "Unmatched:", len(ids_unmatched) #, ", ".join(sorted(ids_unmatched))
    return geneid_to_localization


if __name__ == "__main__":
    main()

