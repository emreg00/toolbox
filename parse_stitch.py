
import gzip, cPickle, os

def main():
    base_dir = "/home/emre/data/drug/stitch/"
    gene_mapping_file = base_dir + "../../interactome/STRING9/entrez_gene_id.vs.string.v9.0.28122012.txt"
    mapping_file = base_dir + "chemicals.inchikeys.v4.0.tsv.gz"
    links_file = base_dir + "9606.protein_chemical.links.detailed.v4.0.tsv.gz"
    output_file = base_dir + "inchi_to_targets.pcl"
    get_inchikey_to_targets(links_file, mapping_file, gene_mapping_file, output_file, cutoff = 900) #, include_score=True)
    return

def get_inchikey_to_targets(links_file, mapping_file, gene_mapping_file, output_file, cutoff = 0, species_prefix = "9606", include_score=False):
    if os.path.exists(output_file):
	inchi_to_targets = cPickle.load(open(output_file))
	return inchi_to_targets 
    # Get inchikey to stitch stereo id mapping
    f = gzip.open(mapping_file)
    line = f.readline()
    cid_to_inchis = {} 
    for line in f:
        #flat_chemical_id stereo_chemical_id source_cid inchikey
	flat_id, stereo_id, source_id, inchikey = line.strip().split()
        #if stereo_id in cid_to_inchi:
        #    print "Overwriting cid:", stereo_id, inchikey, cid_to_inchi[stereo_id]
        #inchi_to_cid[inchikey] = stereo_id
        cid_to_inchis.setdefault(stereo_id, []).append(inchikey)
    f.close()
    # Get ensembl - gene id mapping
    f = open(gene_mapping_file)
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
    f.close()
    # Get inchikey / stitch id to ensembl / gene id mapping
    inchi_to_targets = {}
    f = gzip.open(links_file)
    line = f.readline()
    for line in f:
        #chemical protein experimental prediction database textmining combined_score
	cid, string_id, exp, pred, db, txt, score = line.strip().split()
        if not string_id.startswith(species_prefix):
	    continue
        if cid not in cid_to_inchis:
	    continue
        if string_id not in id_to_geneid:
            continue
        geneid = id_to_geneid[string_id]
	if int(score) >= cutoff:
            for inchi in cid_to_inchis[cid]:
                if include_score:
                    inchi_to_targets.setdefault(inchi, set()).add((geneid, float(score)/1000))
                else:
                    inchi_to_targets.setdefault(inchi, set()).add(geneid)
    f.close()
    f_out = open(output_file, 'w')
    cPickle.dump(inchi_to_targets, f_out)
    f_out.close()
    return inchi_to_targets

def get_cid_to_drugbank_id_mapping(mapping_file, output_file=None):
    if output_file is not None and os.path.exists(output_file):
	cid_to_drugbank_ids = cPickle.load(open(output_file))
	return cid_to_drugbank_ids 
    # Get stitch stereo id to drugbank id mapping
    f = gzip.open(mapping_file)
    line = f.readline()
    cid_to_drugbank_ids = {} 
    for line in f:
        #chemical    alias   source(s)
        words = line.strip().split("\t")
        stereo_id = words[0]
        alias = words[1] 
        sources = words[2].split()
        flag = False
        if alias.startswith("DB"):
            for source in sources:
                if source == "DrugBank":
                    cid_to_drugbank_ids.setdefault(stereo_id, set()).add(alias)
                    flag = True
                    break
            if flag == False:
                print words 
    f.close()
    if output_file is not None:
        f_out = open(output_file, 'w')
        cPickle.dump(cid_to_drugbank_ids, f_out)
        f_out.close()
    return cid_to_drugbank_ids 


if __name__ == "__main__":
    main()

