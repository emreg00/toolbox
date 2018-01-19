
import parse_ncbi

def main():
    mapping_file = "/home/emre/data/ncbi/geneid_to_symbol.txt"
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    file_name = "/home/emre/data/tissue/hpa_subcellular_location.tsv"
    geneid_to_localization = get_geneid_to_localization(file_name, mapping_file)
    print geneid_to_localization["124912"]
    file_name_membrane = "/home/emre/data/tissue/hpa_predicted_membrane.tsv"
    file_name_secreted = "/home/emre/data/tissue/hpa_predicted_secreted.tsv"
    geneid_to_localization = get_geneid_to_predicted_localization(file_name_membrane, file_name_secreted, mapping_file)
    return


def get_geneid_to_localization(file_name, file_name_mapping):
    """
    Relevant extracellular term: Plasma membrane (GO:0005886)
    Some other terms: Vesicles (GO:0043231);Nuclear bodies (GO:0016604);Nucleoplasm (GO:0005654);Cytosol (GO:0005829);Nucleoli (GO:0005730);Mitochondria (GO:0005739);Cell Junctions (GO:0030054);Endoplasmic reticulum (GO:0005783);Centrosome (GO:0005813);Nucleoli fibrillar center (GO:0001650);Golgi apparatus (GO:0005794);Nuclear speckles (GO:0016607)
    """
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    geneid_to_localization = {}
    f = open(file_name)
    f.readline()
    for line in f:
	words = line.strip("\n").split("\t")
	gene, reliability, go = words[1], words[2], words[-1]
	if reliability == "Uncertain":
	    continue
	if gene in name_to_geneid:
	    geneid = name_to_geneid[gene]
	else:
	    print "Unmatched id", gene
	idx = go.find("(")
	go = go[idx+1:-1]
	geneid_to_localization.setdefault(geneid, set()).add(go)
    f.close()
    return geneid_to_localization


def get_geneid_to_predicted_localization(file_name_membrane, file_name_secreted, file_name_mapping):
    """
    Relevant protein classes: Predicted secreted proteins | Predicted membrane proteins
    """
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    geneid_to_localization = {}
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
	    print "Unmatched id", gene
	for localization in protein_class.split(", "):
	    geneid_to_localization.setdefault(geneid, set()).add(localization)
    f.close()
    f = open(file_name_secreted)
    f.readline()
    for line in f:
	words = line.strip("\n").split("\t")
	gene, protein_class, evidence, reliability = words[0], words[6], words[7], words[11]
	if gene in name_to_geneid:
	    geneid = name_to_geneid[gene]
	else:
	    print "Unmatched id", gene
	for localization in protein_class.split(", "):
	    geneid_to_localization.setdefault(geneid, set()).add(localization)
    f.close()
    return geneid_to_localization


if __name__ == "__main__":
    main()

