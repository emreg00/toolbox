
import TsvReader

def main():
    return

def get_geneid_symbol_mapping(file_name):
    """
    To parse Homo_sapiens.gene_info (trimmed to two colums) file from NCBI 
    Creating the file
    wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
    zcat Homo_sapiens.gene_info.gz | cut -f 2,3,5 > geneid_to_symbol.txt
    Remove first line if need be (but keep header)
    """
    f = open(file_name)
    f.readline()
    geneid_to_name = {} # now contains only the official symbol
    name_to_geneid = {}
    for line in f:
	words = line.strip("\n").split("\t")
	if len(words) == 2:
	    geneid, symbol = words
	else:
	    geneid, symbol, alternatives = words[:3]
	    alternatives = alternatives.split("|")
	geneid = geneid.strip() # strip in case mal formatted input file
	symbol = symbol.strip()
	if geneid == "" or symbol == "":
	    continue
	#geneid_to_names.setdefault(geneid, set()).add(symbol) 
	geneid_to_name[geneid] = symbol
	for symbol in [symbol] + alternatives: # added for synonym parsing
	    if symbol in name_to_geneid: 
		if int(geneid) >= int(name_to_geneid[symbol]):
		    continue
		print "Multiple geneids", name_to_geneid[symbol], geneid, symbol
	    name_to_geneid[symbol] = geneid
    f.close()
    return geneid_to_name, name_to_geneid


def get_unigene_to_geneids(file_name, prefix = "Hs."):
    """
    To parse gene2unigene file from NCBI
    """
    f = open(file_name)
    unigene_to_geneids = {}
    f.readline()
    for line in f:
	geneid, unigene = line.strip().split("\t")
	if not unigene.startswith(prefix):
	    continue
	unigene_to_geneids.setdefault(unigene, set()).add(geneid)
    #for unigene, geneids in unigene_to_geneids.iteritems():
    #	if len(geneids) > 1:
    #	    print unigene, geneids
    return unigene_to_geneids


def get_geneid_to_pubmeds(file_name, tax_id = "9606"):
    """
    To parse gene2pubmed file from NCBI 
    """
    f = open(file_name)
    geneid_to_pubmeds = {}
    f.readline()
    for line in f:
	tax, geneid, pubmed_id = line.strip().split("\t")
	if tax != tax_id:
	    continue
	geneid_to_pubmeds.setdefault(geneid, set()).add(pubmed_id)
    return geneid_to_pubmeds


def get_ensembl_protein_to_id(file_name, id_type="geneid"):
    """
    Get ENSEMBL mapping file from BioMart (exporting with ensembl transcipt / protein / gene / symbol / gene id fields)
    id_type: geneid | genesymbol
    """
    ensembl_to_id = {}
    parser = TsvReader.TsvReader(file_name, delim="\t")
    fields_to_include = ["Protein stable ID"]
    if id_type == "geneid":
	column_name = "EntrezGene ID" #"NCBI gene ID" (in GRch38)
    elif id_type == "genesymbol":
	column_name = "HGNC symbol"
    else:
	raise ValueError("Unknown id type")
    fields_to_include += [column_name]
    column_to_index, id_to_values = parser.read(fields_to_include=fields_to_include)
    for ensembl, values in id_to_values.iteritems():
    	for val in values:
    	    gene = val[column_to_index[column_name.lower()]]
	    ensembl_to_id[ensembl] = gene
    return ensembl_to_id


def get_homology_mapping(file_name, tax_id, from_tax_id="9606", symbol_type="geneid"):
    """
    file_name: Homologene data file
    tax_id: Tax id of the species to the genes of which the mapping 
    will be done (e.g., to mouse genes, from human genes)
    from_tax_id: Tax id from which the mapping will be done (default is human)
    symbol_type: geneid | symbol
    Tax ids for popular organisms
    Homo sapiens: 9606
    Mus musculus: 10090
    Rattus norvegicus: 10116
    Drosophila melanogaster: 7227
    Caenorhabditis elegans: 6239
    Saccharomyces cerevisiae: 4932
    Escherichia coli: 562
    Arabidopsis thaliana: 3702
    """
    # 3  9606  34  ACADM  4557231  NP_000007.1
    if symbol_type == "symbol":
	idx = 3
    elif symbol_type == "geneid":
	idx = 2
    else:
	raise ValueError("Unknown symbol type: %s" % symbol_type)
    # Parse group info
    group_to_taxid_to_geneid = {}
    f = open(file_name)
    for line in f:
	words = line.strip("\n").split("\t")
	group, taxid = words[:2]
	geneid = words[idx]
	d = group_to_taxid_to_geneid.setdefault(group, {})
	d[taxid] = geneid
    f.close()
    # Get geneid mapping
    geneid_to_geneid = {}
    for group, taxid_to_geneid in group_to_taxid_to_geneid.iteritems():
	if tax_id not in taxid_to_geneid or from_tax_id not in taxid_to_geneid:
	    continue
	geneid_to_geneid[taxid_to_geneid[from_tax_id]] = taxid_to_geneid[tax_id]
    return geneid_to_geneid, group_to_taxid_to_geneid 

 
