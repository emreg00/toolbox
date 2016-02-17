
import gzip, cPickle, os, re
import parse_string

def main():
    base_dir = "/home/emre/data/"
    chemicals_file = base_dir + "stitch/chemicals.inchikeys.v4.0.1.tsv.gz"
    alias_file = base_dir + "stitch/chemical.aliases.v4.0.tsv.gz"
    links_file = base_dir + "stitch/9606.protein_chemical.links.detailed.v4.0.tsv.gz"
    gene_mapping_file = base_dir + "string/entrez_gene_id.vs.string.v10.28042015.tsv"
    dump_file = base_dir + "stitch/inchi_to_targets.pcl"
    # Get STITCH data
    pubchem_to_drugbank_ids, pubchem_to_target_to_score = get_stitch_info(chemicals_file, alias_file, inchikey_file, links_file, gene_mapping_file, dump_file, species_prefix = "9606")
    # Get inchikey to gene id mapping
    #inchi_to_targets = get_inchikey_to_targets(pubchem_to_target_to_score, pubchem_to_inchi, cutoff = 900)
    # Get drugbank to gene id mapping
    drugbank_id_to_geneids = get_drugbank_to_targets(pubchem_to_drugbank_ids, pubchem_to_target_to_score, cutoff = 900)
    return


def get_stitch_info(chemicals_file, alias_file, inchikey_file, links_file, gene_mapping_file, dump_file, species_prefix = "9606"):
    if os.path.exists(dump_file):
	values = cPickle.load(open(dump_file))
	return values
    # Get cid to drugbank ids
    pubchem_to_drugbank_ids = get_pubchem_to_drugbank_id_mapping(alias_file)
    # Get stitch id to string (ensembl) id mapping
    pubchem_to_target_to_score = get_pubchem_to_targets(links_file, gene_mapping_file, species_prefix)
    f_out = open(dump_file, 'w')
    values = (pubchem_to_drugbank_ids, pubchem_to_target_to_score)
    cPickle.dump(values, f_out)
    f_out.close()
    #pubchem_to_drugbank_ids, pubchem_to_target_to_score # pubchem_to_name_and_smiles, pubchem_to_inchi
    return values 


def get_compound_info_for_given_cids(chemicals_file, inchikey_file, cids):
    # Get cid to name and smiles
    pubchem_to_name_and_smiles = get_pubchem_to_name_and_smiles(chemicals_file, cids) # mem problems
    # Get inchikey to stitch stereo id mapping
    pubchem_to_inchi = get_pubchem_to_inchikeys(inchikey_file, cids) # mem problems
    return pubchem_to_name_and_smiles, pubchem_to_inchi


def get_drugbank_to_targets(pubchem_to_drugbank_ids, pubchem_to_target_to_score, cutoff):
    drugbank_id_to_geneids = {}
    for cid, drugbank_ids in pubchem_to_drugbank_ids.iteritems():
	if cid in pubchem_to_target_to_score:
	    for target, score in pubchem_to_target_to_score[cid].iteritems():
		if score >= cutoff:
		    for drugbank_id in drugbank_ids:
			drugbank_id_to_geneids.setdefault(drugbank_id, set()).add(target)
    return drugbank_id_to_geneids 


def get_pubchem_to_targets(links_file, gene_mapping_file, species_prefix):
    # Get string (ensembl) - gene id mapping
    id_to_geneid = parse_string.get_string_id_to_geneid(gene_mapping_file, species_prefix)
    f = gzip.open(links_file)
    line = f.readline()
    pubchem_to_target_to_score = {}
    for line in f:
        #chemical protein experimental prediction database textmining combined_score
	cid, string_id, exp, pred, db, txt, score = line.strip().split()
        if not string_id.startswith(species_prefix):
	    continue
	if string_id not in id_to_geneid:
	    continue
	cid = cid[3:]
	if cid.startswith("1"):
	    cid = "%s" % (abs(int(cid)) - 100000000)
	else:
	    cid = "%s" % abs(int(cid)) 
	d = pubchem_to_target_to_score.setdefault(cid, {})
	geneid = id_to_geneid[string_id]
	if string_id in d:
	    print "Overwriting", cid, string_id, score
	d[geneid] = int(score) #float(score)/1000
    f.close()
    return pubchem_to_target_to_score


def get_pubchem_to_drugbank_id_mapping(alias_file): #, dump_file=None):
    """
    Parse chemical.aliases file
    """
    #if dump_file is not None and os.path.exists(dump_file):
    #	pubchem_to_drugbank_ids = cPickle.load(open(dump_file))
    #	return pubchem_to_drugbank_ids 
    # Get stitch stereo id to drugbank id mapping
    f = gzip.open(alias_file)
    line = f.readline()
    pubchem_to_drugbank_ids = {} 
    exp = re.compile("DB\d{5}$")
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
                    pubchem_to_drugbank_ids.setdefault(cid, set()).add(alias)
                    flag = True
                    break
                elif source == "845":
                    pubchem_to_drugbank_ids.setdefault(cid, set()).add(alias)
                    flag = True
                    break
            if flag == False:
		if exp.match(alias):
		    print words 
    f.close()
    #if dump_file is not None:
    #    f_out = open(dump_file, 'w')
    #    cPickle.dump(pubchem_to_drugbank_ids, f_out)
    #    f_out.close()
    return pubchem_to_drugbank_ids 


def get_pubchem_to_name_and_smiles(chemicals_file, cids=None):
    """
    Memory intensive, consider parsing it to a ~sql database
    Alternatively, provide a set of cids to be fetched
    """
    f = gzip.open(chemicals_file)
    line = f.readline()
    pubchem_to_name_and_smiles = {}
    for line in f:
	#chemical    name    molecular_weight	SMILES_string
        cid, name, weight, smiles = line.strip().split("\t")
	cid = cid[3:]
	if cid.startswith("1"):
	    cid = "%s" % (abs(int(cid)) - 100000000)
	else:
	    cid = "%s" % abs(int(cid))
	if cids is not None and cid in cids:
	    pubchem_to_name_and_smiles[cid] = (name, smiles)
    f.close()
    return pubchem_to_name_and_smiles


def get_inchikey_to_targets(pubchem_to_target_to_score, pubchem_to_inchi, cutoff):
    inchi_to_targets = {}
    for cid, targets in pubchem_to_target_to_score.iteritems():
        if cid not in pubchem_to_inchi:
	    continue
	for geneid, score in target_to_score.iteritems():
	    if score >= cutoff:
		for inchi in pubchem_to_inchi[cid]:
		    inchi_to_targets.setdefault(inchi, set()).add(geneid)
    return inchi_to_targets


def get_pubchem_to_inchikeys(inchikey_file, cids=None):
    """
    Parse chemicals.inchikeys file
    Memory intensive, consider parsing it to a ~sql database
    Alternatively, provide a set of cids to be fetched
    """
    f = gzip.open(inchikey_file)
    line = f.readline()
    pubchem_to_inchi = {} 
    for line in f:
        #flat_chemical_id stereo_chemical_id source_cid inchikey
	flat_id, stereo_id, source_id, inchikey = line.strip().split()
	cid_flat = "%s" % (abs(int(flat_id[3:])) - 100000000)
	cid_specific = "%s" % abs(int(stereo_id[3:])) 
	if cids is not None:
	    if cid_flat not in cids and cid_specific not in cids:
		continue
        if stereo_id in pubchem_to_inchi:
            print "Overwriting cid:", stereo_id, inchikey, pubchem_to_inchi[stereo_id]
        pubchem_to_inchi.setdefault(cid_flat, set()).add(inchikey)
        pubchem_to_inchi.setdefault(cid_specific, set()).add(inchikey)
    f.close()
    return pubchem_to_inchi


if __name__ == "__main__":
    main()

