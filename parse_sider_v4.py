import TsvReader
import os, cPickle, gzip

def main():
    base_dir = "/home/emre/data/sider"
    indication_file = base_dir + "meddra_all_label_indications.tsv.gz" 
    side_effect_file = base_dir + "meddra_all_label_se.tsv.gz" 
    dump_file = base_dir + "sider.pcl"
    pubchem_to_indications, pubchem_to_side_effects = get_sider_info(indication_file, side_effect_file, dump_file)
    return


def get_sider_info(indication_file, side_effect_file, dump_file):
    if os.path.exists(dump_file):
	pubchem_to_indications, pubchem_to_side_effects = cPickle.load(open(dump_file))
    else:
	pubchem_to_umls_cids, pubchem_to_side_effects = parse_side_effects(side_effect_file)
	pubchem_to_umls_cids, pubchem_to_indications = parse_indications(indication_file)
	cPickle.dump((pubchem_to_indications, pubchem_to_side_effects), open(dump_file, 'w')) 
    return pubchem_to_indications, pubchem_to_side_effects 


def parse_side_effects(side_effect_file):
    f = gzip.open(side_effect_file)
    pubchem_to_umls_cids = {}
    pubchem_to_side_effects = {}
    #cid_to_descriptions = {}
    for line in f.readlines():
	words = line.strip("\n").split("\t")
	# Get Pubchem id (Strip "CID")
	cid_flat = "%s" % (abs(int(words[1][3:])) - 100000000)
	cid_specific = "%s" % abs(int(words[2][3:])) 
	try:
	    term_type, term_cid, term_name = words[4:]
	except:
	    print line
	    print words
	    term_type, term_cid, term_name = words[4:]
	if term_type != "PT":
	    continue
	pubchem_to_umls_cids.setdefault(cid_flat, set()).add(term_cid) 
	pubchem_to_umls_cids.setdefault(cid_specific, set()).add(term_cid) 
	#cid_to_descriptions.setdefault(term_cid, set()).add(term_name.lower())
	pubchem_to_side_effects.setdefault(cid_flat, set()).add(term_name)
	pubchem_to_side_effects.setdefault(cid_specific, set()).add(term_name) 
    f.close()
    return pubchem_to_umls_cids, pubchem_to_side_effects


def parse_indications(indication_file):
    f = gzip.open(indication_file)
    pubchem_to_umls_cids = {}
    pubchem_to_indications = {}
    for line in f.readlines():
	words = line.strip("\n").split("\t")
	# Get Pubchem id (Strip "CID")
	cid_flat = "%s" % (abs(int(words[1][3:])) - 100000000)
	cid_specific = "%s" % abs(int(words[2][3:])) 
	term_type, term_cid, term_name = words[6:]
	if term_type != "PT":
	    continue
	pubchem_to_umls_cids.setdefault(cid_flat, set()).add(term_cid) 
	pubchem_to_umls_cids.setdefault(cid_specific, set()).add(term_cid) 
	pubchem_to_indications.setdefault(cid_flat, set()).add(term_name)
	pubchem_to_indications.setdefault(cid_specific, set()).add(term_name) 
    f.close()
    return pubchem_to_umls_cids, pubchem_to_indications



if __name__ == "__main__":
    main()

