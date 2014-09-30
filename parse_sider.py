import TsvReader
import os, cPickle

def main():
    base_dir = "/home/emre/arastirma/data/drug/sider"
    label_file = base_dir + "label_mapping.tsv"
    side_effect_file = base_dir + "adverse_effects_raw.tsv" 
    dump_file = base_dir + "sider.pcl"
    pubchem_to_siderids, pubchem_specific_to_siderids, siderid_to_descriptions = get_sider_info(label_file, side_effect_file, dump_file)
    return


def get_sider_info(label_file, side_effect_file, dump_file):
    if os.path.exists(dump_file):
	pubchem_to_siderids, pubchem_specific_to_siderids, siderid_to_descriptions = cPickle.load(open(dump_file))
	return pubchem_to_siderids, pubchem_specific_to_siderids, siderid_to_descriptions
    parser = TsvReader.TsvReader(label_file, delim="\t", inner_delim=None)
    column_to_index, id_to_values = parser.read(fields_to_include=["Label", "Name", "Alternative name", "Marker", "CID flat", "CID specific", "URL"], keys_to_include=None, merge_inner_values=False)
    label_to_cids = {}
    label_to_specific_cids = {}
    for label, values in id_to_values.iteritems():
	for val in values:
	    label_to_cids.setdefault(label, set()).add(val[column_to_index["cid flat"]])
	    label_to_specific_cids.setdefault(label, set()).add(val[column_to_index["cid specific"]])
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
    pubchem_to_siderids = {}
    for label, cids in label_to_cids.iteritems():
	for cid in cids:
	    if cid == "":
		continue
	    cid = "%s" % (abs(int(cid)) - 100000000)
	    siderids = pubchem_to_siderids.setdefault(cid, set())
	    if label in label_to_siderids:
		siderids |= label_to_siderids[label]
    pubchem_specific_to_siderids = {}
    for label, cids in label_to_specific_cids.iteritems():
	for cid in cids:
	    if cid == "":
		continue
	    cid = "%s" % abs(int(cid))
	    siderids = pubchem_specific_to_siderids.setdefault(cid, set())
	    if label in label_to_siderids:
		siderids |= label_to_siderids[label]
    cPickle.dump((pubchem_to_siderids, pubchem_specific_to_siderids, siderid_to_descriptions), open(dump_file, 'w'))
    return pubchem_to_siderids, pubchem_specific_to_siderids, siderid_to_descriptions


if __name__ == "__main__":
    main()

