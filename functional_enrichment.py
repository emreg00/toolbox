from funcassociate import client

def main():
    file_name = "gene_list.txt"
    out_file_name = "gene_list_functions.txt"
    f=open(file_name)
    lines = f.readlines()
    a=map(lambda x: x.strip(), lines)
    check_functional_enrichment(a, None, "genesymbol", open(out_file_name, 'w').write, tex_format = False) 
    return


def check_functional_enrichment(subset_gene_ids, gene_ids, id_type, output_method, specie = "Homo sapiens", mode = "unordered", request_info=False, tex_format=False):
    """
	Check GO functional enrichment using funcassociate web service
	gene_ids is a list of gene symbols (without whitespace) or gene ids
	id_type
    """
    if id_type == "geneid":
	id_type = "entrezgene"
    elif id_type == "genesymbol":
	id_type = "hgnc_symbol"
    elif id_type == "uniprotaccession":
	id_type = "uniprot_accession"
    elif id_type == "uniprotentry":
	id_type = "uniprot_id"
    elif id_type == "sgd":
	id_type = "sgd_systematic"
    else:
	raise ValueError("Unrecognized id_type: %s" % id_type)

    reps = 1500
    client_funcassociate = client.FuncassociateClient()
    response = client_funcassociate.functionate(query = subset_gene_ids,
                             species = specie,
                             namespace = id_type,
                             genespace = gene_ids,
                             mode = mode,
                             reps = reps)

    if output_method is None:
	return response["over"]

    #headers = ["N", "M", "X", "LOD", "P", "P_adj", "attrib ID", "attrib name"]
    headers = [ "# of high scoring genes", "# of total genes in high scoring subquery genes", "# of total genes", "Log of odds ratio", "P-value", "Adjusted p-value", "GO term ID", "Go term name" ]

    #if mode == "unordered":
    #	headers.pop(1)
    headers.pop(1) # Now that column is always present independent of the mode
    if tex_format:
	output_method("%s\\\\\n" % " & ".join(headers))
    else:
	output_method("%s\n" % "\t".join(headers))

    zero = "< %f" % (1.0/float(reps))

    for row in response["over"]:
	if mode == "unordered":
	    row = row[:1] + row[2:] #row.pop(1)
        if row[4] is 0:
            row[4] = zero
	if mode == "unordered":
	    interval = range(2,5)
	else:
	    interval = range(3,6)
	#print row
	for i in interval:
	    if isinstance(row[i], str) and row[i].startswith("<"):
		#print row[i]
		val = float(row[i].lstrip("<"))
		if tex_format:
		    row[i] = "$<$%.5f" % val
		else:
		    row[i] = "<%.5f" % val
	    else:
		row[i] = "%.5f" % row[i]
	if tex_format:
	    output_method("%s\\\\\n" % " & ".join(map(str, row)))
	else:
	    output_method("%s\n" % "\t".join(map(str, row)))

    if request_info:
	output_method("\nREQUEST INFO\n")
	info = response["request_info"]
	for k in info.keys():
	    output_method("%s: %s\n" % (k, info[k]))
    return response["over"]


if __name__ == "__main__":
    main()

