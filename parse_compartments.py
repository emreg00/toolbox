
import parse_ncbi

def main():
    file_name_mapping = "/home/emre/data/ncbi/ensembl_id_mapping_grch37.txt"
    ensembl_to_id = parse_ncbi.get_ensembl_protein_to_id(file_name_mapping)
    print ensembl_to_id["ENSP00000269053"]
    file_name = "/home/emre/data/tissue/human_compartment_benchmark.tsv"
    geneid_to_localization = get_geneid_to_localization(file_name, file_name_mapping)
    print geneid_to_localization["124912"]
    return


def get_geneid_to_localization(file_name, file_name_mapping):
    ensembl_to_id = parse_ncbi.get_ensembl_protein_to_id(file_name_mapping)
    geneid_to_localization = {}
    for line in open(file_name):
	ensembl, compartment, sign = line.strip("\n").split("\t")
	if sign == "+":
	    if ensembl in ensembl_to_id:
		geneid = ensembl_to_id[ensembl]
	    #else:
	    #	print "Unmatched id", ensembl
	    geneid_to_localization.setdefault(geneid, set()).add(compartment)
    return geneid_to_localization


if __name__ == "__main__":
    main()

