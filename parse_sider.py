import TsvReader
import os, cPickle

def main():
    base_dir = "/home/emre/arastirma/data/drug/sider"
    label_file = base_dir + "label_mapping.tsv"
    side_effect_file = base_dir + "adverse_effects_raw.tsv" 
    meddra_side_effect_file = base_dir + "meddra_adverse_effects.tsv" 
    dump_file = base_dir + "sider.pcl"
    pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names = get_sider_info(label_file, side_effect_file, dump_file)
    return


def get_sider_info(label_file, side_effect_file, meddra_side_effect_file, dump_file):
    if os.path.exists(dump_file):
	pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names = cPickle.load(open(dump_file))
	return pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names 
    pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names = parse_meddra_side_effects(meddra_side_effect_file)
    #pubchem_to_siderids2, siderid_to_descriptions2, pubchem_to_names2 = parse_side_effects(side_effect_file, label_file)
    cPickle.dump((pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names), open(dump_file, 'w')) 
    return pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names


def parse_meddra_side_effects(meddra_side_effect_file):
    parser = TsvReader.TsvReader(meddra_side_effect_file, delim="\t", inner_delim=None)
    column_to_index, id_to_values = parser.read(fields_to_include=["Name", "CID flat", "CID specific", "SIDER id", "Side effect"], keys_to_include=None, merge_inner_values=False)
    pubchem_to_names = {}
    pubchem_to_siderids = {}
    siderid_to_descriptions = {} # only storing SIDER description
    for name, values in id_to_values.iteritems():
	name = name.lower()
	for val in values:
	    cid_flat, cid_specific, siderid, description = val
	    cid_flat = "%s" % (abs(int(cid_flat)) - 100000000)
	    cid_specific = "%s" % abs(int(cid_specific))
	    pubchem_to_siderids.setdefault(cid_flat, set()).add(siderid) 
	    pubchem_to_siderids.setdefault(cid_specific, set()).add(siderid) 
	    siderid_to_descriptions.setdefault(siderid, set()).add(description.lower())
	    pubchem_to_names.setdefault(cid_flat, set()).add(name)
	    pubchem_to_names.setdefault(cid_specific, set()).add(name) 
    return pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names


def parse_side_effects(side_effect_file, label_file):
    label_to_cids, label_to_specific_cids, label_to_names = parse_labels(label_file)
    label_to_siderids, siderid_to_descriptions = parse_raw_side_effects(side_effect_file)
    pubchem_to_siderids = {}
    pubchem_to_names = {}
    for label, cids in label_to_cids.iteritems():
	for cid in cids:
	    if cid == "":
		continue
	    cid = "%s" % (abs(int(cid)) - 100000000)
	    siderids = pubchem_to_siderids.setdefault(cid, set())
	    if label in label_to_siderids:
		siderids |= label_to_siderids[label]
	    if label in label_to_names:
		for name in label_to_names[label]:
		    pubchem_to_names.setdefault(cid, set()).add(name)
    #pubchem_specific_to_siderids = {}
    for label, cids in label_to_specific_cids.iteritems():
	for cid in cids:
	    if cid == "":
		continue
	    cid = "%s" % abs(int(cid))
	    siderids = pubchem_to_siderids.setdefault(cid, set())
	    #siderids = pubchem_specific_to_siderids.setdefault(cid, set())
	    if label in label_to_siderids:
		siderids |= label_to_siderids[label]
	    if label in label_to_names:
		for name in label_to_names[label]:
		    pubchem_to_names.setdefault(cid, set()).add(name)
    return pubchem_to_siderids, siderid_to_descriptions, pubchem_to_names


def parse_raw_side_effects(side_effect_file):
    parser = TsvReader.TsvReader(side_effect_file, delim="\t", inner_delim=None)
    column_to_index, id_to_values = parser.read(fields_to_include=["Label", "SIDER id", "Side effect"], keys_to_include=None, merge_inner_values=False)
    label_to_siderids = {}
    siderid_to_descriptions = {}
    for label, values in id_to_values.iteritems():
	words = label.split(",")
	for word in words:
	    for val in values:
		siderid, description = val
		label_to_siderids.setdefault(word, set()).add(siderid) 
		siderid_to_descriptions.setdefault(siderid, set()).add(description.lower())
    return label_to_siderids, siderid_to_descriptions


def parse_labels(label_file):
    parser = TsvReader.TsvReader(label_file, delim="\t", inner_delim=None)
    column_to_index, id_to_values = parser.read(fields_to_include=["Label", "Name", "Alternative name", "Marker", "CID flat", "CID specific", "URL"], keys_to_include=None, merge_inner_values=False)
    label_to_cids = {}
    label_to_specific_cids = {}
    label_to_names = {}
    for label, values in id_to_values.iteritems():
	for val in values:
	    label_to_cids.setdefault(label, set()).add(val[column_to_index["cid flat"]])
	    label_to_specific_cids.setdefault(label, set()).add(val[column_to_index["cid specific"]])
	    label_to_names.setdefault(label, set()).add(val[column_to_index["name"]])
    return label_to_cids, label_to_specific_cids, label_to_names


if __name__ == "__main__":
    main()

