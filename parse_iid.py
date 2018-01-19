
import TsvReader, parse_ncbi

def main():
    mapping_file = "/home/emre/data/ncbi/geneid_to_symbol.txt"
    file_name = "/home/emre/data/ppi/iid/iid.human.2017-04.txt"
    tissue = "kidney" #"liver" # None #
    out_file = "/home/emre/data/ppi/iid/network_%s.sif" % tissue
    create_network_file(file_name, mapping_file, out_file, tissue)
    return


def create_network_file(file_name, mapping_file, out_file, tissue=None):
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(mapping_file)
    parser = TsvReader.TsvReader(file_name, delim="\t")
    fields_to_include = ["evidence type", "symbol1", "symbol2"]
    if tissue is not None:
	fields_to_include += [tissue]
    column_to_index, id_to_values = parser.read(fields_to_include=fields_to_include)
    edges = set()
    for evidence, values in id_to_values.iteritems():
	if evidence == "pred":
	    continue
	if tissue is not None:
	    for gene1, gene2, abundance in values:
		#if "P" not in set(abundance):
		if abundance == "-":
		    continue
		if gene1 in name_to_geneid and gene2 in name_to_geneid:
		    geneid1 = name_to_geneid[gene1]
		    geneid2 = name_to_geneid[gene2]
		    edges.add((geneid1, geneid2))
	else:
	    for gene1, gene2 in values:
		if gene1 in name_to_geneid and gene2 in name_to_geneid:
		    geneid1 = name_to_geneid[gene1]
		    geneid2 = name_to_geneid[gene2]
		    edges.add((geneid1, geneid2))
    f = open(out_file, 'w')
    for u, v in edges:
	f.write("%s 1 %s\n" % (u, v))
    f.close()
    return


if __name__ == "__main__":
    main()

