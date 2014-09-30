
def main():
    return

def get_geneid_symbol_mapping(file_name):
    """
    To parse Homo_sapiens.gene_info (trimmed to two colums) file from NCBI 
    """
    f = open(file_name)
    f.readline()
    geneid_to_names = {}
    name_to_geneid = {}
    for line in f:
	geneid, symbol = line.split("\t")
	geneid = geneid.strip()
	symbol = symbol.strip()
	if geneid == "" or symbol == "":
	    continue
	geneid_to_names.setdefault(geneid, set()).add(symbol) 
	if symbol in name_to_geneid: 
	    if int(geneid) >= int(name_to_geneid[symbol]):
		continue
	    print name_to_geneid[symbol], geneid, symbol
	name_to_geneid[symbol] = geneid
    f.close()
    return geneid_to_names, name_to_geneid


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


 
