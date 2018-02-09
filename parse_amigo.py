
import parse_ncbi, parse_uniprot

def main():
    mapping_file = "/home/emre/data/uniprot/idmapping.tab"
    uniprot_to_geneid = parse_uniprot.get_uniprot_to_geneid(mapping_file, uniprot_ids=None)
    file_name = "/home/emre/data/tissue/amigo_extracellular.tsv"
    #file_name = "/home/emre/data/tissue/amigo_membrane.tsv"
    geneid_to_localization = get_geneid_to_localization(file_name, mapping_file)
    print geneid_to_localization["5594"] # MAPK1
    return


def get_geneid_to_localization(file_name, mapping_file):
    name_to_geneid = parse_uniprot.get_uniprot_to_geneid(mapping_file, uniprot_ids=None)
    geneid_to_localization = {}
    uniprots_unmatched = set()
    for line in open(file_name):
	words = line.strip("\n").split("\t")
	uniprot, go, evidence = words[0], words[3], words[7]
	idx = uniprot.find(":")
	uniprot = uniprot[idx+1:]
	#if evidence not in ("EXP", "IDA", "IPA", "IMP", "IGI", "IEP", "HDA"): 
	#    continue
	if uniprot in name_to_geneid:
	    geneid = name_to_geneid[uniprot]
	else:
	    #print "Unmatched id", uniprot
	    if not uniprot.startswith("URS"):
		uniprots_unmatched.add(uniprot)
	    continue
	# Store evidence type as well
	geneid_to_localization.setdefault(geneid, set()).add((go, evidence))
    print "Unmatched:", len(uniprots_unmatched) #, ", ".join(sorted(uniprots_unmatched))
    return geneid_to_localization


def get_geneid_to_extracellular_localization(file_name_membrane, file_name_secreted, mapping_file):
    geneid_to_localization = get_geneid_to_localization(file_name_membrane, mapping_file)
    geneid_to_localization2 = get_geneid_to_localization(file_name_extracellular, mapping_file)
    for geneid, localization in geneid_to_localization2.iteritems():
	geneid_to_localization.setdefault(geneid, set()).add(localization)
    return geneid_to_localization


if __name__ == "__main__":
    main()

