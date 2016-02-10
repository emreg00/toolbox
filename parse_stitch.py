
import gzip, cPickle, os
import parse_string

def main():
    base_dir = "/home/emre/data/"
    chemicals_file = base_dir + "stitch/chemicals.inchikeys.v4.0.1.tsv.gz"
    alias_file = base_dir + "stitch/chemical.aliases.v4.0.tsv.gz"
    links_file = base_dir + "stitch/9606.protein_chemical.links.detailed.v4.0.tsv.gz"
    gene_mapping_file = base_dir + "string/entrez_gene_id.vs.string.v10.28042015.tsv"
    dump_file = base_dir + "stitch/inchi_to_targets.pcl"
    # Get STITCH data
    cid_to_name_and_smiles, cid_to_drugbank_ids, cid_to_target_to_score = get_stitch_info(chemicals_file, alias_file, links_file, dump_file, species_prefix = "9606")
    # Get inchikey to gene id mapping
    #inchi_to_targets = get_inchikey_to_targets(cid_to_target_to_score, cid_to_inchis, cutoff = 900)
    # Get drugbank to gene id mapping
    drugbank_id_to_geneids = get_drugbank_to_targets(cid_to_drugbank_ids, cid_to_target_to_score, cutoff = 900)
    return


def get_stitch_info(chemicals_file, alias_file, links_file, dump_file, species_prefix = "9606"):
    if os.path.exists(dump_file):
	values = cPickle.load(open(dump_file))
	return values
    # Get cid to name and smiles
    cid_to_name_and_smiles = get_cid_to_name_and_smiles(chemicals_file)
    # Get cid to drugbank ids
    cid_to_drugbank_ids = get_cid_to_drugbank_id_mapping(alias_file)
    # Get stitch id to string (ensembl) id mapping
    cid_to_target_to_score = get_cid_to_targets(links_file, species_prefix)
    # Get inchikey to stitch stereo id mapping
    #cid_to_inchis = get_cid_to_inchikeys(inchikey_file)
    f_out = open(dump_file, 'w')
    values = (cid_to_name_and_smiles, cid_to_drugbank_ids, cid_to_target_to_score)
    cPickle.dump(values, f_out)
    f_out.close()
    #cid_to_name_and_smiles, cid_to_drugbank_ids, cid_to_target_to_score
    return values 


def get_drugbank_to_targets(cid_to_drugbank_ids, cid_to_target_to_score, cutoff):
    drugbank_id_to_geneids = {}
    for cid, drugbank_ids in cid_to_drugbank_ids.iteritems():
	if cid in cid_to_target_to_score:
	    for target, score in cid_to_target_to_score[cid].iteritems():
		if score >= cutoff:
		    for drugbank_id in drugbank_ids:
			drugbank_id_to_geneids.setdefault(drugbank_id, set()).add(target)
    return drugbank_id_to_geneids 


def get_cid_to_targets(links_file, gene_mapping_file, species_prefix):
    # Get string (ensembl) - gene id mapping
    id_to_geneid = parse_string.get_string_id_to_geneid(gene_mapping_file, species_prefix)
    f = gzip.open(links_file)
    line = f.readline()
    cid_to_target_to_score = {}
    for line in f:
        #chemical protein experimental prediction database textmining combined_score
	cid, string_id, exp, pred, db, txt, score = line.strip().split()
        if not string_id.startswith(species_prefix):
	    continue
	if string_id not in id_to_geneid:
	    continue
	d = cid_to_target_to_score.setdefault(cid, {})
	d[id_to_geneid[string_id]] = int(score) #float(score)/1000
    f.close()
    return cid_to_target_to_score


def get_cid_to_drugbank_id_mapping(alias_file): #, dump_file=None):
    """
    Parse chemical.aliases file
    """
    #if dump_file is not None and os.path.exists(dump_file):
    #	cid_to_drugbank_ids = cPickle.load(open(dump_file))
    #	return cid_to_drugbank_ids 
    # Get stitch stereo id to drugbank id mapping
    f = gzip.open(alias_file)
    line = f.readline()
    cid_to_drugbank_ids = {} 
    for line in f:
        #chemical    alias   source(s)
        words = line.strip().split("\t")
	cid = words[0][3:]
	if cid.startswith("1"):
	    cid = "%s" % (abs(int(cid)) - 100000000)
	else:
	    cid = "%s" % abs(int(cid)) 
        alias = words[1] 
        sources = words[2].split()
        flag = False
        if alias.startswith("DB"):
            for source in sources:
                if source == "DrugBank":
                    cid_to_drugbank_ids.setdefault(cid, set()).add(alias)
                    flag = True
                    break
            if flag == False:
                print words 
    f.close()
    #if dump_file is not None:
    #    f_out = open(dump_file, 'w')
    #    cPickle.dump(cid_to_drugbank_ids, f_out)
    #    f_out.close()
    return cid_to_drugbank_ids 


def get_cid_to_name_and_smiles(chemicals_file):
    f = gzip.open(chemicals_file)
    line = f.readline()
    cid_to_name_and_smiles = {}
    for line in f:
	#chemical    name    molecular_weight	SMILES_string
        cid, name, weight, smiles = line.strip().split("\t")
	cid = cid[3:]
	if cid.startswith("1"):
	    cid = "%s" % (abs(int(cid)) - 100000000)
	else:
	    cid = "%s" % abs(int(cid)) 
	cid_to_name_and_smiles[cid] = (name, smiles)
    f.close()
    return cid_to_name_and_smiles


def get_inchikey_to_targets(cid_to_target_to_score, cid_to_inchis, cutoff):
    inchi_to_targets = {}
    for cid, targets in cid_to_target_to_score.iteritems():
        if cid not in cid_to_inchis:
	    continue
	for geneid, score in target_to_score.iteritems():
	    if score >= cutoff:
		for inchi in cid_to_inchis[cid]:
		    inchi_to_targets.setdefault(inchi, set()).add(geneid)
    return inchi_to_targets


def get_cid_to_inchikeys(inchikey_file):
    """
    Parse chemicals.inchikeys file
    """
    f = gzip.open(inchikey_file)
    line = f.readline()
    cid_to_inchis = {} 
    for line in f:
        #flat_chemical_id stereo_chemical_id source_cid inchikey
	flat_id, stereo_id, source_id, inchikey = line.strip().split()
	cid_flat = "%s" % (abs(int(flat_id[3:])) - 100000000)
	cid_specific = "%s" % abs(int(stereo_id[3:])) 
        if stereo_id in cid_to_inchi:
            print "Overwriting cid:", stereo_id, inchikey, cid_to_inchi[stereo_id]
        cid_to_inchis.setdefault(cid_flat, set()).add(inchikey)
        cid_to_inchis.setdefault(cid_specific, set()).add(inchikey)
    f.close()
    return cid_to_inchis


if __name__ == "__main__":
    main()

