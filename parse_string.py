
import gzip

def main():
    base_dir = "/home/emre/data/interactome/STRING9/"
    mapping_file = base_dir + "entrez_gene_id.vs.string.v9.0.28122012.txt"
    links_file = base_dir + "protein.links.v9.1.txt.gz"
    output_file = base_dir + "network_500.txt"
    get_interactions(links_file, mapping_file, output_file, cutoff = 500)
    return

def get_interactions(links_file, mapping_file, output_file, cutoff = 0, species_prefix = "9606"):
    f = open(mapping_file)
    line = f.readline()
    id_to_geneid = {} 
    for line in f:
	geneid, stringid = line.strip().split()
	if not stringid.startswith(species_prefix):
	    continue
	if stringid in id_to_geneid: # in case of multiple geneid matches choose lowest
	    if int(geneid) < int(id_to_geneid[stringid]):
		id_to_geneid[stringid] = geneid
	else:
	    id_to_geneid[stringid] = geneid
	#id_to_geneid.setdefault(stringid, set()).add(geneid)
	#if len(id_to_geneid[stringid]) > 1:
	#    print stringid, id_to_geneid[stringid]
    f.close()
    f = gzip.open(links_file)
    f_out = open(output_file, 'w')
    line = f.readline()
    for line in f:
	id1, id2, score = line.strip().split()
	if not id1.startswith(species_prefix) or not id2.startswith(species_prefix):
	    continue
	if id1 not in id_to_geneid or id2 not in id_to_geneid:
	    continue
	if int(score) >= cutoff:
	    f_out.write("%s\t%s\n" % (id_to_geneid[id1], id_to_geneid[id2]))
    f.close()
    f_out.close()
    return 


if __name__ == "__main__":
    main()

