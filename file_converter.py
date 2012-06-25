import network_utilities
import TsvReader

def main():
    input_file = "test.txt"
    mapping_file = "node_mapping.txt.genesymbol"
    output_file = "test.txt.genesymbol"
    #convert_ids_using_mapping_file(input_file, mapping_file, output_file)
    return

def convert_ids_using_mapping_file(input_file, mapping_file, output_file, one_gene_per_node=True, delim="\t"):
    nodes, dummy, node_to_data, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = input_file, store_edge_type = False)

    id_to_mapped_ids = get_id_to_mapped_id_mapping(mapping_file, delim=delim)

    values = []
    #for node, d in node_to_data.iteritems():
    for node in nodes: 
	if node not in id_to_mapped_ids:
	    continue
	if one_gene_per_node:
	    genes = [ id_to_mapped_ids[node][0] ]
	else:
	    genes = id_to_mapped_ids[node]
	for gene in genes:
	    #values.append((d, gene))
	    values.append(gene)
	    
    values.sort()
    values.reverse()
    #i = 1
    f = open(output_file, 'w')
    #f2 = open(output_file + ".ranks", 'w')
    #for d, gene in values:
    for gene in values:
	f.write("%s\n" % (gene))
	#f.write("%s\t%s\n" % (gene, str(score)))
	#f2.write("%s\t%d\n" % (gene, i))
	#i += 1
    f.close()
    #f2.close()
    return

def convert_mapping_file_to_reversed_mapping_file(node_mapping_file, delim="\t"):
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file, delim=delim)
    #f = open(node_mapping_file + ".single", 'w')
    #f.write("user entity id\tgene symbol\n")
    f2 = open(node_mapping_file + ".reversed", 'w')
    words = open(node_mapping_file).readline().strip().split(delim)
    words.reverse()
    f2.write("%s\n" % delim.join(words))
    for node, mapped_ids in id_to_mapped_ids.iteritems():
	#f.write("%s\t%s\n" % (node, mapped_ids[0]))
	f2.write("%s\t%s\n" % (mapped_ids[0], node))
    #f.close()
    f2.close()
    return


def get_id_to_mapped_id_mapping(node_mapping_file, delim="\t", inner_delim = ","):
    reader = TsvReader.TsvReader(node_mapping_file, delim=delim, inner_delim = inner_delim)
    columns, id_to_mapped_ids = reader.read(fields_to_include = None, merge_inner_values = True)
   
    id_to_mapped_ids_formatted = {}
    for node, vals in id_to_mapped_ids.iteritems():
	vals = reduce(lambda x,y: x+y, vals)
	if "-" in vals:
	    vals.remove("-")
	if len(vals) < 1:
	    continue
	id_to_mapped_ids_formatted[node] = vals
    return id_to_mapped_ids_formatted


def create_id_mapping_file_from_gene_info(gene_info_file, gene_ids, output_file):
    reader = TsvReader.TsvReader(gene_info_file, delim="\t", inner_delim = ",")
    columns, id_to_mapped_ids = reader.read(fields_to_include = None, keys_to_include=gene_ids, merge_inner_values = True)

    f = open(output_file, 'w')
    f.write("geneid\tgenesymbol\n")
    for node, vals in id_to_mapped_ids.iteritems():
	vals = reduce(lambda x,y: x+y, vals)
	if "-" in vals:
	    vals.remove("-")
	if len(vals) < 1:
	    continue
	f.write("%s\t%s\n" % (node, vals[0]))
    f.close()
    return

def output_mapped_node_id_scores(output_scores_file, node_mapping_file, one_gene_per_node=True, output_file=None):
    """
	Output mapped ids of nodes 
    """
    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = output_scores_file, store_edge_type = False)

    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)

    values = []
    for node, score in node_to_score.iteritems():
	if node not in id_to_mapped_ids:
	    continue
	if one_gene_per_node:
	    genes = [ id_to_mapped_ids[node][0] ]
	else:
	    genes = id_to_mapped_ids[node]
	for gene in genes:
	    values.append((score, gene))
	    
    values.sort()
    values.reverse()
    included = set()
    i = 1
    if output_file is not None:
	f = open(output_file, 'w')
	f2 = open(output_file + ".ranks", 'w')
	f3 = open(output_file + ".unique", 'w')
	for score, gene in values:
	    f.write("%s\t%s\n" % (gene, str(score)))
	    f2.write("%s\t%d\n" % (gene, i))
	    if gene not in included:
		f3.write("%s\t%s\n" % (gene, str(score)))
	    included.add(gene)
	    i += 1
	f.close()
	f2.close()
	f3.close()
    else:
	print "%s\t%f" % (gene, score)
    return 


if __name__ == "__main__":
    main()


